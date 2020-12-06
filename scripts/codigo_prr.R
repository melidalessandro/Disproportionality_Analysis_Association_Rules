library(tidyverse)

drug <- read.csv("faers_ascii_2020Q1/ascii/DRUG20Q1.txt", sep = "$")
reac <- read.csv("faers_ascii_2020Q1/ascii/REAC20Q1.txt", sep = "$")

drug_reac <- left_join(drug, reac, by = c("primaryid", "caseid"))

drug_reac <- drug_reac %>% 
  filter(prod_ai != "" & prod_ai != "COSMETICS") %>% 
  group_by(prod_ai, pt) %>% 
  tally()

head(drug_reac %>% arrange(desc(n)))

q3 <- quantile(drug_reac$n, 0.75)
q1 <- quantile(drug_reac$n, 0.25)
iqr <- (q3 - q1)
lim_outlier <- (q3 + (1.5 * iqr))

drug_reac_filtered <- drug_reac %>% 
  filter(n >= lim_outlier)


calc_ppr <- function(product, condition) {
  
  tabla <- drug_reac %>% 
    ungroup(prod_ai) %>% 
    mutate(pt = ifelse(pt == condition, "condition", "other_condition"),
           prod_ai = ifelse(prod_ai == product, "product", "other_product"))
  
  tabla <- tabla %>% 
    group_by(pt, prod_ai) %>% 
    summarise(n = sum(n)) %>% 
    pivot_wider(., names_from = pt, values_from = n) %>% 
    arrange(desc(prod_ai))
  
  tabla <- tabla %>% 
    mutate(total_prod = rowSums(tabla[,2:3], na.rm=TRUE))
  
  valor_ppr <- (tabla[1,2] / tabla[1,4]) / (tabla[2,2] / tabla[2,4])
  
  lista <- tibble(product = product, condition = condition, ppr = valor_ppr[1,1])
  
  return(lista)
  
}

library(tictoc)

lista_resultados <- tibble()
drug_reac_filtered <- drug_reac_filtered[1,]

for (i in 1:nrow(drug_reac_filtered)) {
  tic()
  product = drug_reac_filtered[i, 1] %>% pull() %>% as.character()
  condition = drug_reac_filtered[i, 2] %>% pull() %>% as.character()
  tmp_df <- calc_ppr(product, condition)
  lista_resultados <- bind_rows(lista_resultados, tmp_df)
  toc()
}
