library(tidyverse)


#DEMOyyQq.TXT contains patient demographic and administrativeinformation, a single record for each event report.

#DRUGyyQq.TXT contains drug/biologic information for as manymedications as were reported for the event (1 or more per event).

#REACyyQq.TXT contains all "Medical Dictionary for RegulatoryActivities" (MedDRA) terms coded for the event (1 or more). For more

#  information on MedDRA, please contact: TRW, VAR 1/6A/MSSO, 12011 SunsetHills Road, Reston, VA 20190-3285, USA; website is www.meddramsso.com

#OUTCyyQq.TXT contains patient outcomes for the event (0 or more).

#RPSRyyQq.TXT contains report sources for event (0 or more).

#THERyyQq.TXT contains drug therapy start dates and end dates for thereported drugs (0 or more per drug per event).

#INDIyyQq.TXT contains all "Medical Dictionary for RegulatoryActivities" (MedDRA) terms coded for the indications for use(diagnoses) for the reported drugs (0 or more per drug per event).


drug <- read.csv("input/faers_ascii_2020Q1/ASCII/DRUG20Q1.txt", sep = "$")
#demo <- read.csv("input/faers_ascii_2020Q1/ASCII/DEMO20Q1.txt", sep = "$")
#indi <- read.csv("input/faers_ascii_2020Q1/ASCII/INDI20Q1.txt", sep = "$")
#outc <- read.csv("input/faers_ascii_2020Q1/ASCII/OUTC20Q1.txt", sep = "$")
reac <- read.csv("input/faers_ascii_2020Q1/ASCII/REAC20Q1.txt", sep = "$")
#ther <- read.csv("input/faers_ascii_2020Q1/ASCII/THER20Q1.txt", sep = "$")
#rpsr <- read.csv("input/faers_ascii_2020Q1/ASCII/RPSR20Q1.txt", sep = "$")

# Filtro los datos de BI
#datos <- demo %>% 
#  filter(grepl(pattern = "BOEHRINGER INGELHEIM", x = mfr_sndr))

#datos <- left_join(datos, drugs, by = c("primaryid", "caseid")) %>% 
#  left_join(., outc, by = c("primaryid", "caseid")) %>% 
#  left_join(., reac, by = c("primaryid", "caseid"))

#datos_pradaxa1 <- datos %>% 
#  filter(grepl("DABIGATRAN ETEXILATE MESYLATE", prod_ai))
#
#hola <- datos_pradaxa1 %>% distinct(caseid, .keep_all = TRUE)
#
#ggplot(hola, aes(x = sex, y = age)) + 
#  geom_violin() +
#  geom_boxplot(width=0.1)
#
#summary(hola$pt)

drug_reac <- left_join(drug, reac, by = c("primaryid", "caseid"))

drug_count <- drug_reac %>% 
  filter(prod_ai != "" | prod_ai == "COSMETICS") %>% 
  group_by(prod_ai, pt) %>% 
  tally()

ggplot(drug_count) +
  
