# uvozimo podatke
podatki <- read.csv("podatki_1.txt", header = FALSE)

  
# a) Kompletna zadostna statistika
T_1 = sum(log(podatki$V1))
T_2 = sum(log(podatki$V1)^2)


# b) Cenilka največjega verjetja
a_hat = mean(log(podatki$V1))
b2_hat = mean(log(podatki$V1)^2) - mean(log(podatki$V1))^2


# c) Metoda momentov
mu_1 = mean(podatki$V1)
mu_2 = mean(podatki$V1^2)

a_mm = 2*log(mu_1) - 1/2*log(mu_2)
b2_mm = log(mu_2) - 2*log(mu_1)

# d) Metoda delta
a1 = 1/(2*mu_1*mu_2)*(8*mu_1^2 - 20*mu_1*mu_2 + 3*mu_1 + 8*mu_2^2 - 12*mu_2)
a2 = -1/(mu_1*mu_2)*(4*mu_1^2 - 6*mu_1*mu_2 + 2*mu_1 + 2*mu_2^2 - 5*mu_2)
a3 = -2/(mu_1*mu_2)*(2*mu_2 - mu_1)

delta <- matrix(c(a1,a2, a2,a3), nrow=2, ncol=2, byrow=TRUE)
