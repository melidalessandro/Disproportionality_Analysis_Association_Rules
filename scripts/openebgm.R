library(openEBGM)
library(tidyverse)
library(janitor)

options(scipen=999)

#data(caers_raw)
#head(caers_raw, 4)

#drug <- read.csv("input/faers_ascii_2020Q1/ascii/DRUG20Q1.txt", sep = "$", stringsAsFactors = FALSE)
#reac <- read.csv("input/faers_ascii_2020Q1/ascii/REAC20Q1.txt", sep = "$", stringsAsFactors = FALSE)
#drug_reac <- left_join(drug, reac, by = c("primaryid", "caseid"))
#head(drug_reac)

caers_raw <- read_csv("input/CAERS_ASCII_11_18_to_3_20.csv")
caers_raw <- clean_names(caers_raw)
#colnames(caers_raw)

dat <- tidyr::separate_rows(caers_raw, medra_preferred_terms, sep = ", ")

#summary(as.factor(dat$age_units))

# Se pasan las edades a años
dat <- dat %>% 
  mutate(strat_age = case_when(age_units == "Decade(s)" ~ patient_age * 10,
                                 age_units == "day(s)" ~ patient_age / 365,
                                 age_units == "month(s)" ~ patient_age / 12,
                                 age_units == "week(s)" ~ patient_age / 52,
                                 TRUE ~ patient_age),
         strat_age = ifelse(is.na(strat_age), "unknown", ifelse(strat_age < 18, "under_18", "18_plus")),
         sex = ifelse(is.na(sex), "unknown", sex))

# Se cambin los nombres de variables y se seleccionan las que queremos
dat <- dat %>% 
  rename(id = report_id,
         var1 = product,
         var2 = medra_preferred_terms,
         strat_gender = sex) %>% 
  select(id, var1, var2, strat_gender, strat_age) %>% 
  filter(!is.na(var1) & !is.na(var2))

processed <- processRaw(data = dat, stratify = FALSE, zeroes = FALSE)

#resultados <- resultados %>% 
#  mutate(strat_EGBM = ifelse())
#
ggplot(data = processed) +
  geom_point(mapping = aes(x = as.factor(N), y = RR)) +
  scale_y_log10() +
  labs(x = "Cantidad de casos por combinación Producto-Síntoma",
       y = "Log(RR)",
       title = "Varianza de log(RR)")


# Se ahce data squashing para que la hiperoptimizacion sea más rápida

squashed <- squashData(processed) #Using defaults
  
nrow(squashed) / nrow(processed)
#the squashed data set has 10.68% of the observations as the full dataset

squash3 <- autoSquash(processed)
# Esto hace un solo squash por cada N del dataset


#theta_init <- c(alpha1 = 0.2, beta1 = 0.1, alpha2 = 2, beta2 = 4, p = 1/3)
#stats::nlm(negLLsquash, p = theta_init,
#           ni = squashed$N, ei = squashed$E, wi = squashed$weight, N_star = 1)
#
theta_init <- data.frame(alpha1 = c(0.2, 0.1, 0.3),
                         beta1  = c(0.1, 0.2, 0.5),
                         alpha2 = c(2,   10,  6),
                         beta2  = c(4,   10,  6),
                         p      = c(1/3, 0.2, 0.5)
)

squashed2 <- squashData(squashed, count = 2, bin_size = 20, keep_pts = 100)
system.time(
  hyper_estimates_squashed <- autoHyper(data = squashed2, theta_init = theta_init)
)
hyper_estimates_squashed


hyperEM_ests <- hyperEM(squashed, theta_init_vec = c(1, 1, 2, 2, .3),
                        conf_ints = TRUE, track = TRUE)

pdat <- gather(hyperEM_ests$tracking, key = "metric", value = "value", logL:P)
pdat$metric <- factor(pdat$metric, levels = unique(pdat$metric), ordered = TRUE)
ggplot(pdat, aes(x = iter, y = value)) +
  geom_line(size = 1.1, col = "blue") +
  facet_grid(metric ~ ., scales = "free") +
  ggtitle("Convergence Assessment",
          subtitle = "Dashed red line indicates accelerated estimate") +
  labs(x = "Iteration Count", y = "Estimate") +
  geom_vline(xintercept = c(100, 200), size = 1, linetype = 2, col = "red")


#squashed <- squashData(processed)
#squashed2 <- squashData(squashed, count = 2, bin_size = 10, keep_pts = 50)
#theta_init <- data.frame(alpha1 = c(0.2, 0.1, 0.3, 0.5, 0.2),
#                         beta1  = c(0.1, 0.1, 0.5, 0.3, 0.2),
#                         alpha2 = c(2,   10,  6,   12,  5),
#                         beta2  = c(4,   10,  6,   12,  5),
#                         p      = c(1/3, 0.2, 0.5, 0.8, 0.4)
#)
#hyper_estimates <- autoHyper(squashed2, theta_init = theta_init)
#(theta_hat <- hyper_estimates$estimates)

theta_hat <- hyperEM_ests$estimates

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

suspicious <- processed[processed$QUANT_05 >= 2, ]
nrow(suspicious); nrow(processed); nrow(suspicious)/nrow(processed)

suspicious <- suspicious[order(suspicious$QUANT_05, decreasing = TRUE),
                         c("var1", "var2", "N", "E", "QUANT_05", "ebgm", 
                           "QUANT_95")]
head(suspicious, 5)

tabbed <- table(suspicious$var1)
head(tabbed[order(tabbed, decreasing = TRUE)])

#The output above suggests some products which may require further investigation.

ebout <- ebScores(proc, hyperEM_ests = hyper_estimate,
                  quantiles = c(5, 95)) #For the 5th and 95th percentiles
ebout_noquant <- ebScores(proc, hyper_estimate = hyper_estimate,
                          quantiles = NULL) #For no quantiles
