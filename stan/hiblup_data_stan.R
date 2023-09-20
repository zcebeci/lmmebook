# STAN ÖRNEK
# Veri: HIBLUP programındaki testdata.zip dosyaları
# Soyağacı ve fenotip dosyasında IND önekleri kaldırıldı
#
# Paketleri yükle
library(rstan)
library(nadiv)
library(bayesplot)

# Fenotip dosyasını oku
setwd("D:/lmmebook/hiblup/testdata")
pheno <- read.table("phenotype.txt", header=TRUE)
# Soyağacı dosyasını oku
ped <- read.table("pedigree.txt", header=TRUE)

# Akrabalık matrisi
A <- as.matrix(nadiv::makeA(ped))

# X matrisi
# İlk seviyeler olmadan X matrisi
# X <- model.matrix(~sex + season + location + as.numeric(day) + as.numeric(bornweight), data=pheno)
# X <- X[,-1]
#

# X matrisi
# Tüm etki seviyelerini ekleme 
sex <- model.matrix(~-1+sex, data=pheno)
season <- model.matrix(~-1+season, data=pheno)
location <- model.matrix(~-1+location, data=pheno)
bw <- model.matrix(~-1+bornweight, data=pheno)
day <- model.matrix(~-1+day, data=pheno)
X <- cbind(sex, season, location, bw, day)
##################################################

nfix <- ncol(X)
n <- nrow(pheno)
nped <- nrow(ped)

# Z matrisi
Z <- matrix(0, n, nped)
for(i in 1:n)
   Z[i, pheno[i, 1]] <- 1

# Fenotipler
y1 <- as.vector(scale(pheno$tr1))
y2 <- as.vector(scale(pheno$tr2))
y3 <- as.vector(scale(pheno$tr3))

# Stan'e gönderilecek veri
data_stan <- list(J=nfix, N=n, K=nped, X=X, Z=Z, Y=y3, A=A) 

# Stan'i NUTS ile çalıştır

#fit <- stan(file='D:/lmmebook/nishio-arakawa/test.stan', data=data_stan, seed=1234, chain = 1, iter = 100, warmup = 50)
 fit <- stan(file='D:/lmmebook/nishio-arakawa/test.stan', data=data_stan, chain = 1, iter = 500, warmup = 250)
summary(fit)
plot(fit, pars=c("b", "sigma_U", "sigma_E")

# Grafikler ve diyagnostik
library(bayesplot)
color_scheme_set("blue") 
rhats <- rhat(fit, pars=c("b", "sigma_U", "sigma_E"))
rhats
mcmc_rhat(rhats) + yaxis_text(hjust = 1)
