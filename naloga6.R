podatki <- read.delim("podatki_6.txt")
colnames(podatki) <- c("weight", "mpg", "foreign")


# b) beta0, beta1, beta2 po MNV
beta = c(0,0,0)

Z <- podatki[, -3]
Z <- as.matrix(cbind(1, Z))
X <- podatki[,3]

# Newton-Rhapson iz vaj
k = 1

while(k > 0.000001){
  #Z_i * beta za vse i
  z_beta = apply(Z, 1, function(y) y %*% beta)
  p = 1-exp(-exp(z_beta))
  
  #matrika drugih odvodov
  H = t(Z) %*% diag(exp(z_beta)*(X/p - X/p^2 * (1-p) * exp(z_beta) - 1)) %*% Z
  
  A = (X/p - 1) * exp(z_beta)
  #vektor prvih odvodov
  grad = t(Z) %*% A  
  
  #posodobimo beto
  k = solve(H) %*% grad
  beta = beta - k
  k = max(abs(k))
  }

beta



# c) Fisherjeva informacijska matrika
E_H = t(Z) %*% diag(- exp(z_beta)^2 * (1-p)/p) %*% Z
n = length(X)

FI = - E_H/n
FI



# d) Standardne napake za beta
var_beta <- solve(FI)/n
sqrt(var_beta[1,1])
sqrt(var_beta[2,2])
sqrt(var_beta[3,3])



# e) Preizkus domneve H0: beta1 = beta2 = 0
# polni model
l_model = sum(X*log(1-exp(-exp(z_beta))) - (1-X)*exp(z_beta))

# pod H0
beta0 <- c(log(-log(1-mean(X))), 0, 0)
z_beta0 <- apply(Z, 1, function(y) y%*% beta0)
l_H = sum(X*log(1-exp(-exp(z_beta0))) - (1-X)*exp(z_beta0))

# izračunamo testno statistiko po metodi razmerja verjetij
lambda <- -2*(l_H - l_model)
lambda
hi2 <- qchisq(0.95, 2)
hi2

lambda > hi2

# H0 zavrnemo
