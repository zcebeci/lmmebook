# STAN ÖRNEK
# Paketleri yükle
library(rstan)
library(nadiv)
library(bayesplot)

# Fenotip dosyasını oku
pheno <- read.table("D:/lmmebook/datasets/sheep3.dat", header=TRUE)
# Soyağacı dosyasını oku
ped <- read.table("D:/lmmebook/datasets/sheepped3.dat", header=TRUE)

# Akrabalık matrisi
A <- as.matrix(nadiv::makeA(ped))

# X matrisi
cins <- model.matrix(~-1+as.factor(cins), data=pheno)
dtip <- model.matrix(~-1+as.factor(dtip), data=pheno)
dyil <- model.matrix(~-1+as.factor(dyil), data=pheno)
X <- cbind(cins, dtip, dyil)
##################################################

nfix <- ncol(X)
n <- nrow(pheno)
nped <- nrow(ped)

# Z matrisi
Z <- matrix(0, n, nped)
for(i in 1:n)
   Z[i, pheno[i, 1]] <- 1

# Fenotipler
y1 <- as.vector(scale(pheno$typa))
y2 <- as.vector(scale(pheno$lcap))
y3 <- as.vector(scale(pheno$vca))

# Stan'e gönderilecek veri
data_stan <- list(J=nfix, N=n, K=nped, X=X, Z=Z, Y=y3, A=A) 

# Stan'i NUTS ile çalıştır

#fit <- stan(file='D:/lmmebook/nishio-arakawa/test.stan', data=data_stan, seed=1234, chain = 4, iter = 1000, warmup = 500)
 fit <- stan(file='D:/lmmebook/nishio-arakawa/test.stan', data=data_stan, chain = 4, iter = 2500, warmup = 1250)
summary(fit)
plot(fit, pars=c("b", "sigma_U", "sigma_E")

# Grafikler ve diyagnostik
library(bayesplot)
color_scheme_set("blue") 
rhats <- rhat(fit, pars=c("b", "sigma_U", "sigma_E"))
rhats
mcmc_rhat(rhats) + yaxis_text(hjust = 1)
