library(ggplot2)
library(dplyr)


# b) Izračun gamma1 in gamma2
n = 25
p0 = 34/100
C1 = 4
C2 = 13
alpha = 0.05

# funkciji za odvod in vsoto odvodov
odvod <- function(k, n, p){
  choose(n, k) * p^(k - 1) * (1 - p)^(n - k-1) * (k - n*p)
}

vsota_odv <- function(m, M, n, p){
  vsota = 0
  for(i in m:M){
    vsota = vsota + odvod(i, n, p)
  }
  return(vsota)
}

# elementi matrik A in B
a1 = pbinom(C1, n, p0) - dbinom(C1, n, p0)
a2 = 1 - pbinom(C2, n, p0)
a3 = dbinom(C1, n, p0)
a4 = dbinom(C2, n, p0)

b1 = vsota_odv(0, C1 - 1, n, p0)
b2 = vsota_odv(C2 + 1, n, n, p0)
b3 = odvod(C1, n, p0)
b4 = odvod(C2, n, p0)

A <- matrix(c(a3, b3, a4, b4), ncol = 2, nrow = 2)

B <- c(alpha - a1 - a2, - b1 - b2)

gamma <- solve(A, B)
gamma


# c) Pokritost in koef zaupanja
Sp <- c(0.001, 0.002, 0.003, 0.015, 0.033, 0.06, 0.086, 0.115, 0.145, 0.174, 0.206, 0.241,
        0.273, 0.311, 0.347, 0.384, 0.423, 0.465, 0.506, 0.548, 0.593, 0.641, 0.691, 0.742, 0.8, 0.869)
Zg <- rev(rep(1, n+1) - Sp)

N = 10000
pokritost <- rep(0,N)
verjetnosti <- seq(0.001,0.999, length=N)

for(i in 1:N){
  for(j in 1:(n+1)){
    if(Sp[j] <= verjetnosti[i] & verjetnosti[i] <= Zg[j]){
      pokritost[i] <- pokritost[i] + dbinom(j-1, n, verjetnosti[i])
    }}}

koef_zaupanja <- min(pokritost)
koef_zaupanja

graf_pokritosti <- ggplot(data.frame(verjetnosti, pokritost), 
               aes(x = verjetnosti, y = pokritost)) + 
  geom_point(size=0.5, shape=16, color = 'darkblue') +
  xlab("Verjetnosti") + ylab("Pokritost") + ggtitle('Graf pokritosti')
graf_pokritosti




# d) Izračun C1, C2, gamma1, gamma2 za n + 10
n1 = n + 10

# Pregledamo vse možne kombinacije C1 in C2
for(C1 in 0:n1){
  for(C2 in (C1+1):n1){
    # A konstruiramo enako kot v b)
    a3 = dbinom(C1, n1, p0)
    a4 = dbinom(C2, n1, p0)
    b3 = odvod(C1, n1, p0)
    b4 = odvod(C2, n1, p0)
    
    A <- matrix(c(a3, b3, a4, b4), ncol = 2, nrow = 2)
    
    # Robni primeri
    if(C1 == 0 & C2 == n1){
      a1 = 0
      a2 = 0
      b1 = 0
      b2 = 0}
    
    if(C1 == 0 & C2 != n1){
      a1 = 0
      a2 = 1 - pbinom(C2, n1, p0)
      b1 = 0
      b2 = vsota_odv(C2 + 1, n1, n1, p0)}
    
    if(C1 != 0 & C2 == n1){
      a1 = pbinom(C1, n1, p0) - dbinom(C1, n1, p0)
      a2 = 0
      b1 = vsota_odv(0, C1 - 1, n1, p0)
      b2 = 0}
    
    else{
      a1 = pbinom(C1, n1, p0) - dbinom(C1, n1, p0)
      a2 = 1 - pbinom(C2, n1, p0)
      b1 = vsota_odv(0, C1 - 1, n1, p0)
      b2 = vsota_odv(C2 + 1, n1, n1, p0)}
    
    # Matrika B
    B <- c(alpha - a1 - a2, - b1 - b2)
    
    gamma <- tryCatch(solve(A,B), error = function(e) {c(2,2)})

    # preverimo, da gamma1 in gamma2 v [0,1]
    if(gamma[1] >= 0 & gamma[2] >= 0 & gamma[1] <= 1 & gamma[2] <= 1){
      izracuni = c(C1, C2, gamma[1], gamma[2])
    }
  }
}

izracuni



# e) Diagram za n+10
N=1000
p_vektor = seq(0.001, 0.999, length=N)

poisci_C1 <- function(n, p, alpha){
  for (C in 0 : n) {
    if(pbinom(C, n, p) > alpha){
      return(C)
    }
  }}

poisci_C2 <- function(n, p, alpha){
  for (C in 0 : n) {
    if(1-pbinom(C, n, p) < alpha){
      return(C)
    }
  }}

C_min <- lapply(p_vektor, poisci_C1, n=n1, alpha = 0.05/2) %>% unlist()
table(C_min)

C_max <- lapply(p_vektor, poisci_C2, n=n1, alpha = 0.05/2) %>% unlist()
table(C_max)

color <- rep(1:5, (N/5))

df_C <- as.data.frame(cbind(C_min, C_max, p_vektor, color))
df_C$color <- factor(df_C$color)

diagram <- ggplot(df_C) + geom_point(aes(x=p_vektor, y=C_min, colour=color)) + 
  geom_point(aes(x=p_vektor, y=C_max, colour=color)) + 
  scale_color_manual(values = c("darkblue", "darkgreen", "darkred", "darkorange", "purple")) + 
  theme(legend.position = "none") + xlab("Vrednosti p") + ylab("Vrednosti C") + 
  ggtitle("Diagram za n=35") +
  scale_x_continuous(breaks = seq(0,1,0.1)) + 
  scale_y_continuous(breaks = seq(0,35,1))
diagram
