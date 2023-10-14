library(dplyr)
library(partitions)

#H_0
p <- c(1/12, 1/12, 1/3, 1/4, 1/4)

m = 5
#n <- c(40, 60, 80, 100)
n = 40


# Vse možne delitve n na 5 števil
kompozicije <- t(as.matrix(compositions(n,m)))
kompozicije <- as.data.frame(kompozicije)



# a) Velikost preizkusa domene na podlagi razmerja verjetij
chi <- qchisq(0.95, m-1)

test <- kompozicije %>% rowwise() %>% mutate(verjetnost = dmultinom(c(V1,V2,V3,V4,V5), prob = p)) %>% 
  mutate(lambda = prod( (n*p/c(V1,V2,V3,V4,V5))^c(V1,V2,V3,V4,V5))) %>% 
  mutate(statistika = -2*log(lambda)) %>% mutate(test = statistika > chi)

velikost = sum(test[test$test == TRUE,]$verjetnost)

rezultat = paste("Za vzorec velikosti", n, "je velikost preizkusa domneve H_0 na podlagi razmerja verjetij enaka", velikost)
print(rezultat)



# b) Preizkus s Pearsonovo statistiko
#chi <- qchisq(0.95, m-1)

test_Pearson <- kompozicije %>% rowwise() %>% mutate(verjetnost = dmultinom(c(V1,V2,V3,V4,V5), prob = p)) %>% 
  mutate(statistika = sum( (c(V1,V2,V3,V4,V5) - n*p)^2 / (n*p) )) %>% mutate(test = statistika > chi)

velikost_Pearson = sum(test_Pearson[test_Pearson$test == TRUE,]$verjetnost)

rezultat2 = paste("Za vzorec velikosti", n, "je velikost preizkusa domneve H_0 na podlagi Pearsonove statistike enaka", velikost_Pearson)
print(rezultat2)
