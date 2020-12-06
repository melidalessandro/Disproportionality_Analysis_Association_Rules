library(tidyverse)
library(arules)
library(arulesViz)
library(tm)

# Dataset completo

#drug_reac <- drug_reac %>% 
#  distinct(caseid, pt, .keep_all = TRUE)
#
#drug_reac <- drug_reac %>% 
#  mutate(drugname = removePunctuation(toupper(drugname)),
#         pt = removePunctuation(toupper(pt)))
#
#write.csv2(drug_reac, "output/drug_reac_punc.csv", row.names = FALSE)

drug_reac <- read.csv2("output/drug_reac_punc.csv", stringsAsFactors = FALSE)
#
suspicious <- read.csv2("output/suspicious2_ai.csv", stringsAsFactors = FALSE)
#
tabbed <- table(suspicious$var1)
head(tabbed[order(tabbed, decreasing = TRUE)],10)

## PRODUCTOS --------------------------

# Filtar producto
sospechoso <- drug_reac %>% 
  mutate(prod_ai = removePunctuation(toupper(prod_ai))) %>% 
  #  filter(prod_ai == "OCTREOTIDE ACETATE")
  #  filter(prod_ai == "PARATHYROID HORMONE") #RECALLED
 #  filter(prod_ai == "MISOPROSTOL")
  #  filter(prod_ai == "MACITENTAN") ##  LISTO bueno hipertensión arterial pulmonar
  # filter(prod_ai == "OMALIZUMAB") ## bueno asma
  #  filter(prod_ai == "ADALIMUMAB") ## INCORRECT DOSE HEMORRAGIA
  #  filter(prod_ai == "TOCILIZUMAB") ## LISTO 
#filter(prod_ai == "CARBIDOPALEVODOPA") ## bueno, para parkinson
#filter(prod_ai == "EPOPROSTENOL") ##
#filter(prod_ai == "INFLIXIMABDYYB") ## 
filter(prod_ai == "ALEMTUZUMAB") ## 
#filter(prod_ai == "SELEXIPAG") ## 



# Se eliminan los duplicados
sospechoso <- sospechoso %>% 
  distinct(caseid, pt, .keep_all = TRUE)

#write.csv(sospechoso, "output/sospechoso2.csv", row.names = FALSE)

#sospechoso <- read.csv("output/sospechoso2.csv", stringsAsFactors = FALSE)

# Se agrupa por caso y se ven los sintomas
sospechoso_agrup <- sospechoso %>% 
  distinct() %>% 
  select(primaryid, pt) %>% 
  mutate(value = 1) %>% 
  pivot_wider(names_from = pt, values_from = value, values_fn = list(value = length))

### Se eliminan los numeros de casos y se guardan en tid
tid <- as.character(sospechoso_agrup[["primaryid"]])
sospechoso_agrup <- sospechoso_agrup[,-1]

### make all columns factors para pasar a transactions
for(i in 1:ncol(sospechoso_agrup)) sospechoso_agrup[[i]] <- as.factor(sospechoso_agrup[[i]])

trans <- as(sospechoso_agrup, "transactions")

### set transactionID
transactionInfo(trans)[["transactionID"]] <- tid

regla_asoc <- apriori(trans, parameter = list(target="rules", support = 0.05, confidence = 0.50, minlen = 2))
arules::inspect(head(sort(regla_asoc, by="lift", decreasing = TRUE), 50))

subrules <- head(regla_asoc, n = 50, by = "lift")

#plot(regla_asoc, method = "paracoord", cex=0.8, reorder=TRUE)
plot(subrules, method = "graph")

plot(subrules, method = "grouped")

## ARULES con todo

##write.csv(drug_reac_agrup, "output/drug_reac_agrup.csv", row.names = FALSE)
#
#drug_reac_agrup <- read.csv("output/drug_reac_agrup.csv", stringsAsFactors = FALSE)
#
#drug_reac_agrup <- drug_reac %>% 
#  distinct() %>% 
#  select(primaryid, pt) %>% 
#  mutate(value = 1)  
#
#drug_reac_agrup1 <- drug_reac_agrup %>% 
#  slice_sample(n = 50000)
#  
#drug_reac_agrup1 <- drug_reac_agrup1 %>%
#  pivot_wider(names_from = pt, values_from = value, values_fn = list(value = length))
#
## remove transaction IDs
#tid <- as.character(drug_reac_agrup1[["primaryid"]])
#drug_reac_agrup1 <- drug_reac_agrup1[,-1]
#
## make all columns factors
#for(i in 1:ncol(drug_reac_agrup1)) drug_reac_agrup1[[i]] <- as.factor(drug_reac_agrup1[[i]])
#
#trans <- as(drug_reac_agrup1, "transactions")
#
## set transactionID
#transactionInfo(trans)[["transactionID"]] <- tid
#
#regla_asoc <- apriori(trans, parameter = list(target="rules", support = 0.01, confidence = 0.20, minlen = 2))
#arules::inspect(head(sort(regla_asoc, by="lift", decreasing = TRUE), 50))
#
#

## Grafico eventos adversos --------------------------

summary(as.factor(sospechoso$pt))

sospechoso_pt <- sospechoso %>% 
  group_by(pt) %>% 
  tally()

ggplot(data = sospechoso)+
  geom_bar(mapping = aes(pt))







