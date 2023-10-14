library(ggplot2)
library(dplyr)


# uvozimo podatke
podatki <- read.csv("podatki_2.txt", header = FALSE)
n = length(podatki$V1)

# a) Kompletna zadostna statistika
T = sum(podatki$V1)



# c) Najboljši preizkus stopnje znač. 0.05
theta0 = 5*n
alpha = 0.05

# Poiščemo največji C za katerega bo verjetnost >= 1-alpha = 0.95
for (C in 400:800) {
  if(1 - ppois(C, theta0) <= alpha){
    return(C)
  }
}
C

gamma <- (alpha - 1 + ppois(C, theta0))/dpois(C, theta0)
gamma

#1 - ppois(C, theta0) + gamma*dpois(C, theta0)

# aproksimacija z normalno porazdelitvijo
C_2 <- qnorm(1-alpha) * sqrt(theta0) + theta0
C_2

gamma_2 <- (alpha - 1 + ppois(round(C_2), theta0))/dpois(round(C_2), theta0)
gamma_2


theta <- seq(0, 1000, 0.1)
vrednosti <- lapply(theta, function(theta) 1-ppois(C, theta) + gamma * dpois(C, theta)) %>% unlist()

graf_moci <- ggplot(data.frame(theta, vrednosti), aes(theta, vrednosti)) + geom_line( color = 'darkblue') +
  xlab("Theta") + ylab("Moč") + ggtitle('Graf funkcije moči') + 
  geom_vline(aes(xintercept = theta0), color ="red", lty = "dashed") +
  geom_hline(aes(yintercept = alpha), color = "black", lty = "dashed") +
  annotate(geom="text", x=50, y=0.08, label="alpha = 0.05", col="black") + 
  annotate(geom="text", x=400, y=1, label="theta0 = 500", col="red") 

graf_moci



# d) Interval zaupanja

for (theta in seq(0, 1000, 0.1)) {
  if(ppois(T, theta) <= 1 - alpha/2){
    theta_inf = theta
    break()
  }}
for (theta in seq(1000, 0, -0.1)) {
  if(ppois(T, theta) >= alpha/2){
    theta_sup = theta
    break()
  }}
  
IZ_ntheta = c(theta_inf, theta_sup)

IZ_theta = IZ_ntheta/100
