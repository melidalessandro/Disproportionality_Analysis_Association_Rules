library(tidyverse)

# Predicción --------------------

#drug_reac <- read.csv2("input/drug_reac.csv", stringsAsFactors = FALSE)
#
#drug_reac <- drug_reac %>% 
#  distinct(caseid, pt, .keep_all = TRUE) %>% 
#  mutate_all(na_if,"")
#
#write.csv(drug_reac, "input/drug_reac_distinct.csv", row.names = FALSE)

# Drug - Reac ------------------------------
# Importo el archivo menos pesado sin duplicados
drug_reac <- read.csv2("input/drug_reac_distinct.csv", stringsAsFactors = FALSE)

#summary(as.factor(drug_reac$cum_dose_unit))
colnames(drug_reac)

data <- drug_reac %>% 
  select(primaryid, drug_seq, role_cod, drugname, prod_ai, route,
         dose_vbm, cum_dose_chr, cum_dose_unit, dechal,
         rechal, exp_dt, dose_amt, dose_unit, dose_form,
         dose_freq, pt)
# drug_rec_act se saca porque tiene muchos NA
## Se deben normalizar algunas variables acá como dosis

# Demograficos -------------------------
demo <- read.delim("input/faers_ascii_2020Q1/ascii/DEMO20Q1.txt", sep = "$", encoding = "UTF-8")

# Preprocesamiento de datos demograficos
demo <- demo %>% 
  mutate_all(na_if,"") %>% 
  mutate(wt = ifelse(wt_cod != "LBS", wt, 0.453592 * wt),
         age = as.numeric(age),
         age = case_when(age_cod == "DEC" ~ age * 10,
                         age_cod == "DY" ~ age / 365,
                         age_cod == "HR" ~ age / (24 * 365),
                         age_cod == "MON" ~ age / 12,
                         age_cod == "WK" ~ age / 52,
                         TRUE ~ age),
         age = ifelse(is.na(age), "UNK", age),
         sex = ifelse(is.na(sex)| sex == "UNK", "UNK", sex)) %>% 
  select(primaryid, event_dt, age, sex, wt, rept_dt,occr_country)

data <- left_join(data, demo, by = "primaryid")

# Outcome ---------------------------
outc <- read.delim("input/faers_ascii_2020Q1/ascii/OUTC20Q1.txt", sep = "$", encoding = "UTF-8")

outc <- outc %>% 
  select(primaryid, outc_cod)

data <- left_join(data, outc, by = "primaryid")
# Se agregan duplicados porque para una misma droga puede haber varios outcomes

colnames(data)

# Therapy ---------------------------
ther <- read.delim("input/faers_ascii_2020Q1/ascii/THER20Q1.txt", sep = "$", encoding = "UTF-8")
# el join se debe hacer por primary id y drug_seq con dsg_drug_seq
# puede haber duplicados porque puede haber empezado y terminado un tratamiento

# Hay que calcular la duracion y agrupar por pirmay id sumando duracion de tratamientos
#ther <- ther %>% 
#  mutate()

# OPSUMIT ------------------------

data1 <- data %>% 
  filter(drugname == "OPSUMIT")

data1 %>%  group_by(primaryid) %>% tally()

summary(as.factor(data1$cum_dose_unit))

data1_results <- data1 %>% 
  select(primaryid, outc_cod)