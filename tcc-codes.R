library(tidyverse)
library(dplyr)
library(lubridate)
library(broom)
library(readr)
library(readxl)
library(fpc)
library(factoextra)
library(cluster)
library(plotly)
library(minpack.lm)
library(NbClust)
library(patchwork) 
library(lmtest)
library(nlme)
library(MASS)

# Leitura dos dados

covid <- read_csv("caso_full.csv", col_types = cols(city_ibge_code = col_character()))

pop22 <- read_csv("pop22.csv", col_types = cols(state_code = col_character())) %>% 
    mutate(city_ibge_code = paste0(state_code,city_ibge_code)) %>%
    dplyr::select(codmun = city_ibge_code, pop2022 = pop)

regiao <- read.csv("estados.csv") %>% dplyr::select(regiao, uf)

area <- read_csv("area.csv", col_types = cols(codmun = col_character())) %>% dplyr::select(codmun, area)

cnes <- read_csv("cnes.csv", col_types = cols(codmun = col_character()))

profissionais <- read_csv("profissionais.csv", col_types = cols(codmun = col_character()))


# Preparação da base

cvd_window <- covid %>% filter(!is.na(city), !is.na(state), !is.na(city_ibge_code)) %>%
  distinct() %>% 
  group_by(city, state, city_ibge_code) %>% 
  dplyr::select(estado = state,municipio =  city, pop = estimated_population, data = date, new_confirmed, codmun = city_ibge_code) %>% 
  mutate(new_confirmed = if_else(new_confirmed < 0, 0, new_confirmed),
         casosAcumulado = cumsum(new_confirmed)) %>% 
  mutate(K = max(casosAcumulado) + 1) %>%
  filter(casosAcumulado > 0) %>%
  filter(data < (min(data) + 30)) %>% 
  filter(n() >= 30) %>%
  mutate(max = max(casosAcumulado), 
    min = min(casosAcumulado),
    Range = max - min, 
    Ratio = max/min) %>%
  ungroup() %>%
  arrange(municipio, estado, codmun, data)  %>%  
  inner_join(pop22, join_by(codmun == codmun)) %>%
  inner_join(area, join_by(codmun == codmun)) %>%
  inner_join(regiao, join_by(estado == uf)) %>%
  mutate(codmun = substring(codmun,1,6)) %>%
  inner_join(cnes, join_by(codmun == codmun)) %>%
  inner_join(profissionais, join_by(codmun == codmun)) %>%
  mutate(cnes_per_pop = 1000*cnes/pop2022) %>%
  mutate(prof_per_pop = 1000*profissionais/pop2022) %>%
  mutate(densidade = pop2022/area)

# Análise Exploratoria

graph = cvd_window %>% group_by(municipio, estado, codmun, regiao, pop2022) %>% 
                       summarise(casosAcumulado = max(casosAcumulado)) %>% ungroup() %>%
                       mutate(perc_casos = 10000*casosAcumulado/pop2022)

summary(graph$perc_casos)

View(graph %>% group_by(regiao) %>% summarise(media = mean(perc_casos), mediana = median(perc_casos)))

qplot(perc_casos, data = graph, geom = "histogram",
      xlab = "Casos de COVID-19 nos primeiros 30 dias a cada 1000 habitantes", ylab = "Frequência",  
      color = I("black"), fill = I("wheat"))

cvd_window %>%
  group_by(municipio, estado, codmun) %>% 
  count() %>% 
  arrange(n) %>%
  print(n=5)

cvd_window %>%
  dplyr::select(municipio, estado, codmun) %>% 
  distinct %>% count()

# Informações básicas

capitais = c(280030, 150140, 310620, 140010, 530010, 500270, 510340, 410690, 420540, 230440,
             520870, 250750, 160030, 270430, 130260, 240810, 172100, 431490, 110020, 261160,
             120040, 330455, 292740, 211130, 355030, 221100, 320530)

basic_infos = cvd_window %>% 
  mutate(capitais = case_when(codmun %in% capitais ~ 1, TRUE ~ 0)) %>% dplyr::select(codmun, densidade, 
                                    cnes_per_pop, prof_per_pop, regiao, municipio, estado, capitais) %>% distinct

# Preparação para a Regressão

cvd_data <- cvd_window %>% 
  dplyr::select(regiao, estado, municipio, codmun, data, casosAcumulado, K) %>% 
  group_by(codmun) %>% 
  mutate(Days = data - min(data) + 1, 
         numDays = as.numeric(Days), 
         exponential = log(casosAcumulado),
         gompertz = log(log(K/casosAcumulado))
  )

# Regressões

regressions <- cvd_data %>% 
  group_by(codmun, municipio, estado) %>% 
  nest() %>%   
  mutate(
    fit_linear = map(data, ~ lm(casosAcumulado ~ numDays, data = .x)),
    glanced_linear = map(fit_linear, glance),
    tidied_linear = map(fit_linear, tidy),
    augmented_linear = map(fit_linear, augment),
    
    fit_log = map(data, ~ lm(exponential ~ numDays, data = .x)),
    glanced_log = map(fit_log, glance),
    tidied_log = map(fit_log, tidy),
    augmented_log = map(fit_log, augment),
    
    fit_loglog = map(data, ~ lm(gompertz ~ numDays, data = .x)),
    glanced_loglog = map(fit_loglog, glance),
    tidied_loglog = map(fit_loglog, tidy),
    augmented_loglog = map(fit_loglog, augment)
  )

# Exemplos

## Linear

max_r.squared <- regressions %>% 
  ungroup() %>% 
  unnest(glanced_log) %>% 
  filter(adj.r.squared == max(adj.r.squared))

model_linear <- regressions %>% 
  ungroup() %>% 
  unnest(tidied_linear) %>% 
  filter(codmun == 210520) %>% 
  dplyr::select(term, estimate) %>% 
  pivot_wider(names_from = term, values_from = estimate)

cvd_data %>% 
  filter(codmun == 210520) %>%
  mutate(fitted = numDays*model_linear$numDays + model_linear$`(Intercept)`,
         modelo = 'Linear')  %>% 
  ggplot() +
  geom_point(aes(x = numDays, y = casosAcumulado)) +
  theme_light() +
  geom_line(aes(x = numDays, y = fitted, color = modelo), size = 1)+
  labs(title = "Igarapé Grande - MA", x = "Número de dias a partir do primeiro caso", y = "Casos acumulados de COVID-19") 

## Exponencial

max_r.squared <- regressions %>% 
  ungroup() %>% 
  unnest(glanced_log) %>% 
  mutate(adj.r.squared = case_when(is.na(adj.r.squared) ~ 0, TRUE ~ adj.r.squared)) %>%
  dplyr::select(municipio,estado,codmun, adj.r.squared)

View(max_r.squared)

max_r.squared <- regressions %>% 
  ungroup() %>% 
  unnest(glanced_log) %>% 
  mutate(adj.r.squared = case_when(is.na(adj.r.squared) ~ 0, TRUE ~ adj.r.squared)) %>%
  filter(adj.r.squared == max(adj.r.squared))

model_log <-regressions %>% 
  ungroup() %>% 
  unnest(tidied_log) %>% 
  filter(codmun == 320530) %>% 
  dplyr::select(term, estimate) %>% 
  pivot_wider(names_from = term, values_from = estimate)


cvd_data %>% 
  filter(codmun == 320530) %>%
  mutate(fitted = exp(numDays*model_log$numDays + model_log$`(Intercept)`), 
         modelo = "Exponencial") %>% 
  ggplot() +
  geom_point(aes(x = numDays, y = casosAcumulado)) +
  theme_light() +
  geom_line(aes(x = numDays, y = fitted, color = modelo), size = 1)+
  labs(title = "Vitória - ES",  x = "Número de dias a partir do primeiro caso", y = "Casos acumulados de COVID-19")

## Gompertz

max_r.squared <- regressions %>% 
  ungroup() %>% 
  unnest(glanced_loglog) %>% 
  mutate(adj.r.squared = case_when(is.na(adj.r.squared) ~ 0, TRUE ~ adj.r.squared)) %>%
  dplyr::select(municipio,estado,codmun, adj.r.squared)

View(max_r.squared)

max_r.squared <- regressions %>% 
  ungroup() %>% 
  unnest(glanced_log) %>% 
  mutate(adj.r.squared = case_when(is.na(adj.r.squared) ~ 0, TRUE ~ adj.r.squared)) %>%
  filter(adj.r.squared == max(adj.r.squared))

model_loglog <-regressions %>% 
  ungroup() %>% 
  unnest(tidied_loglog) %>% 
  filter(codmun == 520870) %>% 
  dplyr::select(term, estimate) %>% 
  pivot_wider(names_from = term, values_from = estimate)


cvd_data %>% 
  filter(codmun == 520870) %>%
  mutate(fitted = K/exp(exp(numDays*model_loglog$numDays + model_loglog$`(Intercept)`)),
         modelo = 'Gompertz') %>% 
  ggplot() +
  geom_point(aes(x = numDays, y = casosAcumulado)) +
  theme_light() +
  geom_line(aes(x = numDays, y = fitted, color = modelo), size = 1)+
  labs(title = "Goiânia - GO",  x = "Número de dias a partir do primeiro caso", y = "Casos acumulados de COVID-19")

## Recife

cvd_data %>% filter(municipio == 'Recife')
max_r.squared <- regressions %>% 
  ungroup() %>% 
  unnest(tidied_loglog) %>% 
  filter(term == "numDays") %>% 
  filter(estimate == min(estimate))

model_linear <-regressions %>% 
  ungroup() %>% 
  unnest(tidied_linear) %>% 
  filter(codmun == 261160) %>% 
  dplyr::select(term, estimate) %>% 
  pivot_wider(names_from = term, values_from = estimate)

model_log <-regressions %>% 
  ungroup() %>% 
  unnest(tidied_log) %>% 
  filter(codmun == 261160) %>% 
  dplyr::select(term, estimate) %>% 
  pivot_wider(names_from = term, values_from = estimate)

model_loglog <-regressions %>% 
  ungroup() %>% 
  unnest(tidied_loglog) %>% 
  filter(codmun == 261160) %>% 
  dplyr::select(term, estimate) %>% 
  pivot_wider(names_from = term, values_from = estimate)

linear <- cvd_data %>% 
  filter(codmun == 261160) %>% 
  mutate(fitted = numDays*model_linear$numDays + model_linear$`(Intercept)`,
         modelo = 'Linear') 

gompert <- cvd_data %>% 
  filter(codmun == 261160) %>% 
  mutate(fitted = K/exp(exp(numDays*model_loglog$numDays + model_loglog$`(Intercept)`)),
         modelo = 'Gompertz')

cvd_data %>% 
  filter(codmun == 261160) %>% 
  mutate(fitted = exp(numDays*model_log$numDays + model_log$`(Intercept)`), 
         modelo = "Exponencial") %>% 
  bind_rows(gompert) %>% 
  bind_rows(linear) %>% 
  ggplot() +
  geom_point(aes(x = numDays, y = casosAcumulado)) +
  theme_light() +
  geom_line(aes(x = numDays, y = fitted, color = modelo), size = 1)+
  labs(title = 'Recife - PE',  x = "Número de dias a partir do primeiro caso", y = "Casos acumulados de COVID-19") 


# Análise de Resíduo

calculate_tests <- function(fit_model, data) {
  residuals <- residuals(fit_model)
  
  if (length(unique(residuals)) == 1) {
    return(data.frame(Shapiro = NA, GQ = NA))
  }
  
  shapiro_p_value <- shapiro.test(residuals)$p.value
  gq_p_value <- gqtest(fit_model, fraction = 0.5)$p.value
  
  return(data.frame(Shapiro = shapiro_p_value, GQ = gq_p_value))
}

# Aplicar a função para cada município
results_log <- regressions %>% 
  rowwise() %>% 
  do(calculate_tests(.$fit_log, .$data))

results_linear <- regressions %>% 
  rowwise() %>% 
  do(calculate_tests(.$fit_linear, .$data))

results_loglog <- regressions %>% 
  rowwise() %>% 
  do(calculate_tests(.$fit_loglog, .$data))

# Histograma dos p-valores de Shapiro-Wilk
par(mfrow = c(1,3))

hist(results_linear$Shapiro, main="Linear", xlab="p valor", ylab = "Frequência", breaks=50)
abline(v = 0.05, col = "red")

hist(results_log$Shapiro, main="Exponencial", xlab="p valor", ylab = "Frequência", breaks=50)
abline(v = 0.05, col = "red")

hist(results_loglog$Shapiro, main="Gompertz", xlab="p valor", ylab = "Frequência", breaks=50)
abline(v = 0.05, col = "red")

# Histograma dos p-valores de Goldfeld-Quandt
par(mfrow = c(1,3))

hist(results_linear$GQ, main="Linear", xlab="p valor", ylab = "Frequência", breaks=50)
abline(v = 0.05, col = "red")

hist(results_log$GQ, main="Exponencial", xlab="p valor", ylab = "Frequência", breaks=50)
abline(v = 0.05, col = "red")

hist(results_loglog$GQ, main="Gompertz", xlab="p valor", ylab = "Frequência", breaks=50)
abline(v = 0.05, col = "red")


# Residuals medio

avg_residuals_linear <- regressions %>%
  unnest(augmented_linear) %>%
  group_by(numDays) %>%
  summarise(mean_residual = mean(.resid))

avg_residuals_log <- regressions %>%
  unnest(augmented_log) %>%
  group_by(numDays) %>%
  summarise(mean_residual = mean(.resid))

avg_residuals_loglog <- regressions %>%
  unnest(augmented_loglog) %>%
  group_by(numDays) %>%
  summarise(mean_residual = mean(.resid))

par(mfrow = c(1,3))
ggplot(avg_residuals_linear, aes(x = numDays, y = mean_residual)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_minimal() +
  labs(
       x = "Número de Dias",
       y = "Média dos Resíduos")

ggplot(avg_residuals_log, aes(x = numDays, y = mean_residual)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_minimal() +
  labs(
    x = "Número de Dias",
    y = "Média dos Resíduos")

ggplot(avg_residuals_loglog, aes(x = numDays, y = mean_residual)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_minimal() +
  labs(
    x = "Número de Dias",
    y = "Média dos Resíduos")

par(mfrow = c(1,3))
qqnorm(avg_residuals_linear$mean_residual, xlab = "Quantis teóricos da distribuição normal", 
       ylab = "Quantis observados da amostra", main = "")
qqline(avg_residuals_linear$mean_residual)

qqnorm(avg_residuals_log$mean_residual, xlab = "Quantis teóricos da distribuição normal", 
       ylab = "Quantis observados da amostra", main = "")
qqline(avg_residuals_log$mean_residual)

qqnorm(avg_residuals_loglog$mean_residual, xlab = "Quantis teóricos da distribuição normal", 
       ylab = "Quantis observados da amostra", main = "")
qqline(avg_residuals_loglog$mean_residual)

# Clustering

vazios = cvd_window %>% filter(Range == 0) %>% dplyr::select(codmun) %>% mutate(vazio = 1) %>% distinct()

clusterization = regressions %>% 
  ungroup() %>% 
  
  unnest(glanced_linear) %>% 
  dplyr::select(rsquared_linear = r.squared, pvalue_linear = p.value, glanced_log, glanced_loglog, codmun) %>% 
  
  unnest(glanced_log) %>%
  dplyr::select(rsquared_linear, pvalue_linear, rsquared_log = r.squared,  pvalue_log = p.value, glanced_loglog, codmun) %>% 
  
  unnest(glanced_loglog) %>%
  dplyr::select(rsquared_linear, pvalue_linear, rsquared_log, pvalue_log, 
         rsquared_loglog = r.squared, pvalue_loglog = p.value, codmun)

## Distribuição dos municipios por p-valor

clusterization_distrib <- clusterization %>%
  left_join(vazios, join_by(codmun == codmun)) %>%
  mutate(
    rsquared_linear = case_when(vazio == 1 ~ 0, TRUE ~ rsquared_linear),
    rsquared_log = case_when(vazio == 1 ~ 0, TRUE ~ rsquared_log),
    rsquared_loglog = case_when(vazio == 1 ~ 0, TRUE ~ rsquared_loglog))

clusterization_distrib  %>% mutate(
  pvalue_linear = case_when(vazio == 1 ~ 1, pvalue_linear < 0.05 ~ 1, pvalue_linear >= 0.05 ~ 0),
  pvalue_log = case_when(vazio == 1 ~ 1, pvalue_log < 0.05 ~ 1, pvalue_log >= 0.05 ~ 0),
  pvalue_loglog = case_when(vazio == 1 ~ 1, pvalue_loglog < 0.05 ~ 1, pvalue_loglog >= 0.05 ~ 0)) %>%
  group_by(pvalue_linear,pvalue_log, pvalue_loglog) %>% count

## Manipulações e ajustes

clusterization_test <- clusterization %>%
  left_join(vazios, join_by(codmun == codmun)) %>%
  mutate(
    rsquared_linear = case_when(vazio == 1 ~ 0, TRUE ~ rsquared_linear),
    rsquared_log = case_when(vazio == 1 ~ 0, TRUE ~ rsquared_log),
    rsquared_loglog = case_when(vazio == 1 ~ 0, TRUE ~ rsquared_loglog))  %>%
  mutate(
    pvalue_linear = case_when(vazio == 1 ~ 1, pvalue_linear < 0.05 ~ 1, pvalue_linear >= 0.05 ~ 0),
    pvalue_log = case_when(vazio == 1 ~ 1, pvalue_log < 0.05 ~ 1, pvalue_log >= 0.05 ~ 0),
    pvalue_loglog = case_when(vazio == 1 ~ 1, pvalue_loglog < 0.05 ~ 1, pvalue_loglog >= 0.05 ~ 0)) %>%
  mutate(
    rsquared_linear = case_when(pvalue_linear < 0.05 ~ 0, TRUE ~ rsquared_linear),
    rsquared_log = case_when(pvalue_log < 0.05 ~ 0, TRUE ~ rsquared_log),
    rsquared_loglog = case_when(pvalue_loglog < 0.05 ~ 0, TRUE ~ rsquared_loglog)) %>%
    dplyr::select(rsquared_linear, rsquared_log, rsquared_loglog, codmun)

clusterization_test$vazio = NULL

clusterization_final = clusterization_test %>% column_to_rownames(var = "codmun")

summary(clusterization_final)

# Metódo para definição dos tamanhos de clusters

fviz_nbclust <- function (x, FUNcluster = NULL, method = c("silhouette", "wss", 
                                                           "gap_stat"), diss = NULL, k.max = 10, nboot = 100, verbose = interactive(), 
                          barfill = "steelblue", barcolor = "steelblue", linecolor = "steelblue", 
                          print.summary = TRUE, ...) 
{
  set.seed(123)
  if (k.max < 2) 
    stop("k.max must bet > = 2")
  method = match.arg(method)
  if (!inherits(x, c("data.frame", "matrix")) & !("Best.nc" %in% 
                                                  names(x))) 
    stop("x should be an object of class matrix/data.frame or ", 
         "an object created by the function NbClust() [NbClust package].")
  if (inherits(x, "list") & "Best.nc" %in% names(x)) {
    best_nc <- x$Best.nc
    if (any(class(best_nc) == "numeric") ) 
      print(best_nc)
    else if (any(class(best_nc) == "matrix") )
      .viz_NbClust(x, print.summary, barfill, barcolor)
  }
  else if (is.null(FUNcluster)) 
    stop("The argument FUNcluster is required. ", "Possible values are kmeans, pam, hcut, clara, ...")
  else if (!is.function(FUNcluster)) {
    stop("The argument FUNcluster should be a function. ", 
         "Check if you're not overriding the specified function name somewhere.")
  }
  else if (method %in% c("silhouette", "wss")) {
    if (is.data.frame(x)) 
      x <- as.matrix(x)
    if (is.null(diss)) 
      diss <- stats::dist(x)
    v <- rep(0, k.max)
    if (method == "silhouette") {
      for (i in 2:k.max) {
        clust <- FUNcluster(x, i, ...)
        v[i] <- .get_ave_sil_width(diss, clust$cluster)
      }
    }
    else if (method == "wss") {
      for (i in 1:k.max) {
        clust <- FUNcluster(x, i, ...)
        v[i] <- .get_withinSS(diss, clust$cluster)
      }
    }
    df <- data.frame(clusters = as.factor(1:k.max), y = v, 
                     stringsAsFactors = TRUE)
    ylab <- "Total Within Sum of Square"
    if (method == "silhouette") 
      ylab <- "Average silhouette width"
    p <- ggpubr::ggline(df, x = "Grupos", y = "y", group = 1, 
                        color = linecolor, ylab = ylab, xlab = "Número de grupos k")
    if (method == "silhouette") 
      p <- p + geom_vline(xintercept = which.max(v), linetype = 2, 
                          color = linecolor)
    return(p)
  }
  else if (method == "gap_stat") {
    extra_args <- list(...)
    gap_stat <- cluster::clusGap(x, FUNcluster, K.max = k.max, 
                                 B = nboot, verbose = verbose, ...)
    if (!is.null(extra_args$maxSE)) 
      maxSE <- extra_args$maxSE
    else maxSE <- list(method = "firstSEmax", SE.factor = 1)
    p <- fviz_gap_stat(gap_stat, linecolor = linecolor, 
                       maxSE = maxSE)
    return(p)
  }
}

.viz_NbClust <- function (x, print.summary = TRUE, barfill = "steelblue", 
                          barcolor = "steelblue") 
{
  best_nc <- x$Best.nc
  if (any(class(best_nc) == "numeric") )
    print(best_nc)
  else if (any(class(best_nc) == "matrix") ) {
    best_nc <- as.data.frame(t(best_nc), stringsAsFactors = TRUE)
    best_nc$Number_clusters <- as.factor(best_nc$Number_clusters)
    if (print.summary) {
      ss <- summary(best_nc$Number_clusters)
      cat("Among all indices: \n===================\n")
      for (i in 1:length(ss)) {
        cat("*", ss[i], "proposed ", names(ss)[i], 
            "as the best number of clusters\n")
      }
      cat("\nConclusion\n=========================\n")
      cat("* According to the majority rule, the best number of clusters is ", 
          names(which.max(ss)), ".\n\n")
    }
    df <- data.frame(Number_clusters = names(ss), freq = ss, 
                     stringsAsFactors = TRUE)
    p <- ggpubr::ggbarplot(df, x = "Number_clusters", 
                           y = "freq", fill = barfill, color = barcolor) + 
      labs(x = "Número de grupos k", y = "Frequência entre todos os índices", 
           title = paste0("O número de grupos ótimo é k = ", 
                          names(which.max(ss))))
    return(p)
  }
}

environment(fviz_nbclust) <- asNamespace("factoextra")
assignInNamespace("fviz_nbclust",fviz_nbclust,"factoextra")
environment(.viz_NbClust) <- asNamespace("factoextra")
assignInNamespace(".viz_NbClust",.viz_NbClust,"factoextra")

# NbClust

nb <- NbClust(clusterization_final, distance = "euclidean", min.nc = 2,
              max.nc = 6, method = "kmeans")

fviz_nbclust(nb)

# K-means

set.seed(16)
km.res <- kmeans(clusterization_final, 3, nstart = 100, iter.max = 100)

# Biplot PCA

fviz_cluster(        km.res,
                     data=clusterization_final,
                     stand = FALSE,
                     show.clust.cent = FALSE,
                     geom = 'point',
                     title = ''
)

# Análise dos clusters

agg = clusterization_test
agg$codmun = NULL
aggregate(agg, by=list(cluster=km.res$cluster), median)
aggregate(agg, by=list(cluster=km.res$cluster), mean)

clusters = regressions %>% ungroup() %>% 
 inner_join(basic_infos, join_by(codmun == codmun)) %>%
 mutate(cluster = km.res$cluster)


clusters %>% group_by(cluster) %>% summarise(mean(densidade), mean(cnes_per_pop), mean(prof_per_pop))
clusters %>% group_by(cluster) %>% summarise(median(densidade), median(cnes_per_pop), median(prof_per_pop))


clusters %>% group_by(cluster) %>% count() %>% print()
clusters %>% group_by(cluster, regiao) %>% count() %>% arrange(cluster, -n) %>% print()
clusters %>% group_by(cluster, capitais) %>% count() %>% arrange(cluster, -n) %>% print()
