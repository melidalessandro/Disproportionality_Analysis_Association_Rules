library(GGally)

demo <- read.delim("input/faers_ascii_2020Q1/ascii/DEMO20Q1.txt", sep = "$", encoding = "UTF-8", stringsAsFactors = FALSE)

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
         sex = ifelse(is.na(sex)| sex == "UNK", "UNK", sex))

sosp <- left_join(sospechoso, demo, by = "primaryid")

sosp1 <- sosp %>% 
  select(drugname, prod_ai, mfr_sndr, age, sex, occr_country) %>% 
  mutate(mfr_sndr = as.factor(mfr_sndr))

hola <- sosp1 %>% 
  select(mfr_sndr) %>% 
  distinct()

summary(as.factor(sosp$pt))

ggpairs(sosp1 %>% select(-occr_country, -prod_ai))
