library(tidyverse)
library(biclust)

#processed <- read.csv2("output/processed2_ai.csv", stringsAsFactors = FALSE)

# Se agrupa por caso y se ven los sintomas
#biclust <- processed %>% 
#  select(var1, var2, ebgm) %>% 
#  mutate(ebgm = ifelse(ebgm > 5, 1, 0)) %>% 
#  pivot_wider(names_from = var2, values_from = ebgm, values_fill = 0) 
#  
#biclust0 <- biclust %>% 
#  mutate_if(is.numeric, ~replace(., is.na(.), 0))
#
#rownames(biclust0) <- biclust0$var1
#
#biclust0 <- biclust0[,-1]
#write.csv2(biclust0, "output/biclust1.csv", row.names = FALSE)

biclust <- read.csv2("output/biclust1.csv")

x <- as.matrix(biclust)

x <- discretize(biclust)
res <- biclust(x, method=BCBimax(), minr=2, minc=2, number=100)
res

parallelCoordinates(x = x, bicResult=res, number=100)

drawHeatmap(x = x, bicResult=res, number = 1)

drawHeatmap2(x = x, bicResult = res, number = 1) # shown in picture below
drawHeatmap2(x = x, bicResult = res, number = 2)
