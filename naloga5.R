podatki <- read.delim("podatki_5.txt")


# funkciji za normo in GS algoritem iz vaj
norma = function(x){
  sqrt(sum(x^2))
}


GS = function(Z){
  n <- nrow(Z)
  d <- ncol(Z)
  S <- matrix(0, nrow = n, ncol = d)
  P <- matrix(0, nrow = d, ncol = d)
  # prvi stolpec S samo normiramo
  # drop = F uporabljamo zato, da se ohrani oblika (rezultat
  # bo še vedno matrika in ne vektor)
  S[, 1] = Z[, 1, drop = F] / norma(Z[, 1])
  P[1, 1] = norma(Z[, 1])
  
  for (i in 2:d){
    S[, i] = Z[, i, drop = F]
    # odštevamo pravokotne projekcije
    for (j in 1:(i-1)){
      koef = sum(Z[, i] * S[, j]) # skalarni produkt
      S[, i] = S[, i] - koef * S[, j]
      P[j, i] = koef
    }
    normaS = norma(S[,i])
    # če je pripadajoč stolpec enak 0, bo tudi norma enaka 0
    if (normaS < 1e-10){ 
      # na diagonalo P lahko postavimo katerokoli pozitivno število,
      # ker se množi z 0 (preko S)
      P[i, i] = 1
    }else{
      S[, i] = S[, i, drop = F] / normaS
      P[i, i] = normaS
    }
  }
  # vrnemo seznam z S in P
  list(S, P)
}


#--------------------------------------------------------------------------------
# a) Predpostavimo model LDL = beta0 + beta1*TCH + beta2*HDL + beta3*TRI + eps

#(i) Ocena beta0, beta1, beta2, beta3 po MNK
Z <- podatki[c(1:3)]
Z <- as.matrix(cbind(Z0 = 1, Z))
X <- as.matrix(podatki[4])

# Posplošeni inverz
S = GS(Z)[[1]]
P = GS(Z)[[2]]

J = t(S) %*% S
J = round(J, digits = 10)

# cenilka s posplošenim inverzom
inv_splosni = solve(P) %*% J %*% solve(t(P))

beta_hat_splosni = inv_splosni %*% t(Z) %*% X
beta_hat_splosni

# cenilka s klasičnim inverzom
beta_hat_klas = solve(t(Z) %*% Z) %*% t(Z) %*% X
beta_hat_klas

# cenilka s funkcijo lm
lm = lm(LDL ~ TCH + HDL + TRI, data=podatki) 
lm$coeff



# (ii) Preizkus domneve H0: beta1 = beta2 = beta3 = 0 pri alpha = 0.05
alpha = 0.05
n = nrow(Z)
d = ncol(Z)
r = 3

# Cenilka za beto pod domnevo H0
beta_hat_H <- c(mean(X), 0, 0, 0)

# Vsota kvadratov residualov v izhodiščnem modelu
VKR <- t(X - Z %*% beta_hat_klas) %*% (X - Z %*% beta_hat_klas) 

# Statistika
F_stat <- as.numeric(((t(Z %*% beta_hat_klas - Z %*% beta_hat_H) %*% (Z %*% beta_hat_klas - Z %*% beta_hat_H)) / r ) /
  (VKR/(n - d)))

#test
fisher = qf(1-alpha, df1 = r, df2 = n - d)
F_stat > fisher 

# Ker velja > hipotezo zavrnemo



# (iii) Preizkus domneve beta0 = 0
#alpha = 0.05

# Statistika
l <- c(1, 0, 0, 0)
t_stat <- as.numeric(abs(t(beta_hat_splosni) %*% l / ( sqrt(t(inv_splosni %*% l) %*% l) %*% sqrt(VKR / (n-d)))))

#test
student = qt(1-alpha/2, df = n - d)
t_stat > student 

# H0 ponovno zavrnemo


# Statistiki pri (ii) in (iii) lahko preberemo tudi iz modela lm
F_stat2 <- as.numeric(summary(lm)$fstatistic[1])

beta_df <- data.frame(summary(lm)$coeff)
t_stat2 <- beta_df['(Intercept)' , 't.value']


#--------------------------------------------------------------------------------
# b) Predpostavimo model LDL = beta1*TCH + beta2*HDL + beta3*TRI + eps

# (i) Ocena beta1, beta2, beta3 po MNK
Z <- as.matrix(podatki[c(1:3)])

# Posplošeni inverz
S = GS(Z)[[1]]
P = GS(Z)[[2]]

J = t(S) %*% S
J = round(J, digits = 10)

# cenilka s posplošenim inverzom
inv_splosni = solve(P) %*% J %*% solve(t(P))

beta_hat_splosni = inv_splosni %*% t(Z) %*% X
beta_hat_splosni

# cenilka s klasičnim inverzom
beta_hat_klas = solve(t(Z) %*% Z) %*% t(Z) %*% X
beta_hat_klas

# cenilka s funkcijo lm
lm = lm(LDL ~0 + TCH + HDL + TRI, data=podatki) 
lm$coeff



# (ii) Preizkus domneve H0: (beta1, beta2, beta3) = (1, -1, -0.45)
# alpha = 0.05
r = 3
d = ncol(Z)
beta_hat_H <- c(1, -1, -0.45)

# Vsota kvadratov residualov v izhodiščnem modelu
VKR <- t(X - Z %*% beta_hat_klas) %*% (X - Z %*% beta_hat_klas) 

# Statistika
F_stat <- as.numeric(((t(Z %*% beta_hat_klas - Z %*% beta_hat_H) %*% (Z %*% beta_hat_klas - Z %*% beta_hat_H)) / r ) /
                       (VKR/(n - d)))

#test
fisher = qf(1-alpha, df1 = r, df2 = n - d)
F_stat > fisher 

# H0 spet zavrnemo



# (iii) Območje zaupanja za (Beta1, Beta2, Beta3), stopnja zaupanja 0.95

# Diagonaliziramo matriko (Z^T)*Z
A <- t(Z) %*% Z

Q <- eigen(A)$vectors
Q

lambda <- diag(eigen(A)$values)
lambda

#sum(abs(Q %*% lambda %*% t(Q) - A)) < 1e-10

# izračunamo še ksi
ksi = fisher * (d / (n - d)) * VKR
ksi

# dolžine polosi
a = as.numeric(ksi/sqrt(lambda[1,1]))
b = as.numeric(ksi/sqrt(lambda[2,2]))
c = as.numeric(ksi/sqrt(lambda[3,3]))

V_elipsoid = 4/3 * pi * a * b * c



# (iv) Bonferronijev popravek
t_bon = qt(1 - 0.05/(2*d), n - d)

# Meje za intervale
delta1 = t_bon * sqrt(inv_splosni[1, 1] * VKR /(n - d))
delta2 = t_bon * sqrt(inv_splosni[2, 2] * VKR /(n - d))
delta3 = t_bon * sqrt(inv_splosni[3, 3] * VKR /(n - d))

# Intervali zaupanja
IZ1 <- c(beta_hat_splosni[1] - delta1, beta_hat_splosni[1] + delta1)
IZ2 <- c(beta_hat_splosni[2] - delta2, beta_hat_splosni[2] + delta2)
IZ3 <- c(beta_hat_splosni[3] - delta3, beta_hat_splosni[3] + delta3)


V_kvader = (IZ1[2] - IZ1[1]) * (IZ2[2] - IZ2[1]) * (IZ3[2] - IZ3[1])
