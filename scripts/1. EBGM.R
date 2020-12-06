library(openEBGM)
library(tidyverse)
library(stringi)
library(stringr)
library(tm)
library(glue)

## Calculo EBGM -----------------------------------------------------------

options(scipen = 9999)

#drug <- read.delim("input/faers_ascii_2020Q1/ascii/DRUG20Q1.txt", sep = "$", encoding = "UTF-8")
#reac <- read.delim("input/faers_ascii_2020Q1/ascii/REAC20Q1.txt", sep = "$", encoding = "UTF-8")
#drug_reac <- left_join(drug, reac, by = c("primaryid", "caseid"))

drug_reac <- read.csv2("input/drug_reac.csv", stringsAsFactors = FALSE)

#drug_reac %>% distinct(pt) %>% tally()

drug_reac <- drug_reac %>% 
  distinct(caseid, pt, .keep_all = TRUE) %>% 
  select(primaryid, drugname, pt)

colnames(drug_reac)

drogadrug_reacs_1 <- drug_reac %>% 
  distinct(prod_ai) #2659

demo <- read.delim("input/faers_ascii_2020Q1/ascii/DEMO20Q1.txt", sep = "$", encoding = "UTF-8")

demo <- demo %>% 
  mutate_all(na_if,"") %>% 
  mutate(wt = ifelse(wt_cod != "LBS", wt, 0.453592 * wt),
         strat_age = as.numeric(age),
         strat_age = case_when(age_cod == "DEC" ~ strat_age * 10,
                               age_cod == "DY" ~ strat_age / 365,
                               age_cod == "HR" ~ strat_age / (24 * 365),
                               age_cod == "MON" ~ strat_age / 12,
                               age_cod == "WK" ~ strat_age / 52,
                               TRUE ~ strat_age),
         strat_age = ifelse(is.na(strat_age), "UNK", ifelse(strat_age < 18, "under_18", "18_plus")),
         sex = ifelse(is.na(sex)| sex == "UNK", "UNK", sex)) %>% 
  select(primaryid, strat_age, sex)

caers_raw <- left_join(drug_reac, demo, by = "primaryid")

drogas <- caers_raw %>% 
  distinct(drugname) #12992

drogas_1 <- caers_raw %>% 
  distinct(prod_ai) #2659

eventos <- caers_raw %>% 
  distinct(pt) #12160

caers_raw <- caers_raw %>% 
  mutate(drugname = removePunctuation(toupper(drugname)),
         pt = removePunctuation(toupper(pt)))

drogas <- caers_raw %>% 
  distinct(drugname) #12129

eventos <- caers_raw %>% 
  distinct(pt) #12160

rm(demo, drug_reac)

# Se cambin los nombres de variables y se seleccionan las que queremos
dat <- caers_raw %>% 
  mutate_all(na_if,"") %>% 
  rename(id = primaryid,
         var1 = prod_ai,
         var2 = pt,
         strat_gender = sex) %>% 
  select(id, var1, var2, strat_gender, strat_age) %>% 
  drop_na()

#summary(as.factor(caers_raw$prod_ai))

processed <- processRaw(data = dat, stratify = FALSE, zeroes = FALSE)

#write.csv2(processed, "output/processed_sinprocesar_ai.csv", row.names = FALSE)

processed <- read.csv2("output/processed_sinprocesar_ai.csv", stringsAsFactors = FALSE)

# Se ahce data squashing para que la hiperoptimizacion sea más rápida

squashed <- squashData(processed) #Using defaults

nrow(squashed) / nrow(processed)
#the squashed data set has 36% of the observations as the full dataset
# ahora 43%

#squash3 <- autoSquash(processed)
# Esto hace un solo squash por cada N del dataset
#nrow(squash3) / nrow(processed)
#colnames(squashed)

theta_init <- c(alpha1 = 0.2, beta1 = 0.1, alpha2 = 2, beta2 = 4, p = 1/3)
hyper_ests <- stats::nlm(negLLsquash, p = theta_init,
                         ni = squashed$N, ei = squashed$E, wi = squashed$weight, N_star = 1)

theta_hat <- hyper_ests$estimate
theta_hat
# 0.000001022186 0.024754918795 2.049054840029 3.963300939862 0.429404051046
#0.00000003693272 0.02514681354248 2.05502850467530 3.95917745274168 0.45198930394566

squashed2 <- squashData(squashed, count = 2, bin_size = 100, keep_pts = 1)
#nrow(squashed2) / nrow(processed)

#squashed3 <- squashData(squashed2, count = 3, bin_size = 100, keep_pts = 1)
#squashed4 <- squashData(squashed3, count = 4, bin_size = 100, keep_pts = 1)
#squashed5 <- squashData(squashed4, count = 5, bin_size = 100, keep_pts = 1)

#system.time(
#  hyper_estimates_squashed <- autoHyper(data = squashed2, theta_init = theta_init,
#                                        zeroes = FALSE, squashed = TRUE, N_star = 1, max_pts = 100000)
#)
#hyper_estimates_squashed
#
#hyperEM_ests <- hyperEM(squash3, theta_init_vec = c(0.1, 0.1, 2, 4, .3),
#                        conf_ints = TRUE, track = TRUE)

## OCon squashed 2

hyper_ests_2 <- stats::nlm(negLLsquash, p = theta_init,
                           ni = squashed2$N, ei = squashed2$E, wi = squashed2$weight, N_star = 1)

theta_hat_2 <- hyper_ests_2$estimate
theta_hat_2
#0.0000002252281 0.0059366418874 3.8037183809090 2.7935776987590 0.1869345942400 #inicial
#0.0000002782572 0.0240538756247 2.0486022185356 3.9637744364176 0.4201357357919 #squashed_2
#0.0000008896992 0.0247379022936 2.0489887929265 3.9633519718162 0.4295803153719
#0.000002449404 0.025138795527 2.054990060185 3.959196892011 0.452170479833
qn <- Qn(theta_hat, N = processed$N, E = processed$E)
head(qn)

identical(length(qn), nrow(processed))

summary(qn)

processed$ebgm <- ebgm(theta_hat, N = processed$N, E = processed$E, qn  = qn)
head(processed)

processed$QUANT_05 <- quantBisect(5, theta_hat = theta_hat,
                                  N = processed$N, E = processed$E, qn = qn)
processed$QUANT_95 <- quantBisect(95, theta_hat = theta_hat,
                                  N = processed$N, E = processed$E, qn = qn)
head(processed)

#write.csv2(processed, "output/processed2_ai.csv", row.names = FALSE)

processed <- read.csv2("output/processed2_ai.csv", stringsAsFactors = FALSE)

summary(processed)

## EBGM > 5 -------------------------------------------

suspicious <- processed[processed$QUANT_05 >= 2, ]
nrow(suspicious); nrow(processed); nrow(suspicious)/nrow(processed)

nrow(suspicious)/nrow(processed)

head(suspicious)

suspicious <- suspicious[order(suspicious$QUANT_05, decreasing = TRUE),
                         c("var1", "var2", "N", "E", "PRR", "QUANT_05", "ebgm", 
                           "QUANT_95")]
#write.csv2(suspicious, "output/suspicious2_ai.csv", row.names = FALSE)

suspicious <- read.csv2("output/suspicious2_ai.csv", stringsAsFactors = FALSE)

suspicious <- suspicious %>% 
  mutate(log_PRR = log10(PRR))

head(suspicious, 5)

tabbed <- table(suspicious$var1)
head(tabbed[order(tabbed, decreasing = TRUE)])

#The output above suggests some products which may require further investigation.

tabbed1 <- table(suspicious$var1) %>% as.data.frame() %>% filter(Freq > 10)
head(tabbed[order(tabbed, decreasing = TRUE)])


## Gráfico varianza log RR-------------------------------------

processed_grafico <- processed %>% 
  mutate(strat_egbm = ifelse(ebgm > 5, TRUE, FALSE),
         log_PRR = log10(PRR)) %>% 
  filter(log_PRR < Inf)

#summary(processed_grafico)

processed_grafico_2 <- processed_grafico %>%
  slice(0:5000) 

#processed_grafico <- read.csv2("output/processed_grafico_ai.csv", stringasFactors = FALSE)

ggplot(data = processed_grafico) +
  geom_point(mapping = aes(x = N, y = log_PRR, color = strat_egbm), alpha = 0.3) +
  labs(x = "Cantidad de casos por combinación Droga-Evento",
       y = "Log(PRR)",
       title = "Variación de log(PRR) por cantidad de combinaciones Droga-Evento",
       color = "EBGM > 5") +
  scale_x_discrete(limits = c(seq(0,100,20)))+
  xlim(0,100) +
  theme_minimal()

#write.csv2(processed_grafico, "output/processed_grafico_ai.csv", row.names = FALSE)

## Gráfico calamar -------------------------------------------

processed_cal <- processed %>% 
  mutate(N_modif = ifelse(N >= 9, "9+", N),
         log_PRR = log(PRR)) %>% 
  filter(log_PRR < Inf)

processed_cal_2 <- processed_cal %>%
  slice(0:5000) 

ggplot(data = processed_cal) +
  geom_point(mapping = aes(x = log_PRR, y = ebgm, color = as.factor(N_modif))) +
  labs(x = "Log(PRR)",
       y = "EBGM",
       title = "Encogimiento de EBGM",
       color = "N") +
  theme_minimal()

## Histograma EBGM ----------------------

ggplot(processed, aes(x = ebgm)) + 
  geom_histogram(binwidth = 0.5, fill = "#F8766D", origin = 0.5) +
  # ylim(0,50000) +
  xlim(0,8) +
  labs(x = "EBGM",
       y = "Frecuencia",
       title = "Histograma de EBGM") +
  theme_minimal()

ggplot(processed %>% filter(ebgm >= 5), aes(x = ebgm)) + 
  geom_histogram(binwidth = 0.5, fill = "#7CAE00", origin = 0.5) +
  # ylim(0,50000) +
  xlim(5,50) +
  labs(x = "EBGM",
       y = "Frecuencia",
       title = "Histograma de EBGM >= 5") +
  theme_minimal()
  
# Como el histograma no queda bien, hago una tabla con frecuencias

histograma <- as.data.frame(table(cut(processed$ebgm, breaks=seq(0,max(processed$ebgm), by=5))))

#write.csv2(histograma, "output/histograma.csv", row.names = FALSE)

## Grafico de intervalos de confianza -------------

drug <- suspicious %>% distinct(var1)

event <- suspicious %>% distinct(var2)

sosp_20 <- suspicious %>% 
  slice(1:20) %>% 
  mutate(var1 = ifelse(var1 == "COAGULATION FACTOR IX HUMANCOAGULATION FACTOR VII HUMANCOAGULATION FACTOR X HUMANPROTEIN CPROTEIN S HUMANPROTHROMBIN", "COAGULATION FACTOR IX (*)", var1),
         droga_evento = glue("{var1} - {var2}"),
         N_modif = glue("N = {N}"))

# Esto se pone para que se ordene el valor de ebgm en el grafico
sosp_20$droga_evento = factor(sosp_20$droga_evento, levels=sosp_20[order(sosp_20$ebgm), "droga_evento"])

ggplot(data = sosp_20, mapping = aes(x = droga_evento, y = ebgm, colour = droga_evento)) +
  geom_point() +
  geom_errorbar(aes(x = droga_evento,
                    ymin = QUANT_05, ymax = QUANT_95,
                    colour = droga_evento),
                size = 0.75,
                width = 0.4,                    # Width of the error bars
                position=position_dodge(0.9),
                inherit.aes = FALSE) +
  geom_text(aes(label = N_modif), hjust = 0.50, vjust = 1.25,
            size = 3.5) +
  coord_flip() +
  labs(title = "Top 20 EBGM para Drogas-Eventos",
       x = "Combinación Droga-Evento",
       y = "Valor de EBGM",
       caption = "(*) COAGULATION FACTOR IX HUMANCOAGULATION FACTOR VII HUMANCOAGULATION FACTOR X HUMANPROTEIN CPROTEIN S HUMANPROTHROMBIN") +
  theme_minimal() +
  theme(legend.position = "none")

suspicious[18,1]

hormona <- processed %>% filter(var1 == "PARATHYROID HORMONE")

## SUSPICIOUS -----------------------

var2_hist <- suspicious %>% 
  group_by(var2) %>% 
  tally() %>% 
  arrange(desc(n)) %>% 
  slice(0:10)

# Esto se pone para que se ordene el valor de ebgm en el grafico
var2_hist$n = factor(var2_hist$n, levels=sort(table(var2_hist$n)))

ggplot(var2_hist, aes(x = reorder(var2, n), y = n, fill = var2)) + 
  geom_col() +
  coord_flip() +
  geom_text(aes(label = n), hjust = 2, vjust = 0.25,
            size = 3.5) +
  labs(x = "Eventos Adversos",
       y = "Frecuencia",
       title = "Top 10 eventos adversos para casos sospechosos") +
  theme_minimal() +
  theme(legend.position = "none")

