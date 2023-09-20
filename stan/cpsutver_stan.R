# Stan programı / Ceylanpınar Tarım İşletmesi 1995 yılı 1. laktasyon süt verimi
#
library(rstan)
library("QTLRel")
library(nadiv)
library(bayesplot)
library(ggplot2)

cpsutver <- read.table("D:/lmmebook/datasets/cpsutver2.dat",header=T)
cpsutver$hayvan <- as.factor(cpsutver$hayvan)
cpsutver$baba <- as.factor(cpsutver$baba)
cpsutver$ana <- as.factor(cpsutver$ana)
cpsutver$bay <- as.factor(cpsutver$bay)
cpsutver$byas <- as.numeric(cpsutver$byas)
cpsutver$sver305 <- as.numeric(cpsutver$sver305)
str(cpsutver)
head(cpsutver, 3)

# Eklemeli genetik akrabalık matrisi
cpped <- cpsutver[,1:3]
colnames(cpped) <- c("id","father","mother")
cpped <- QTLRel::pedRecode(cpped)

n <- nrow(cpsutver)    # Fenotipi olan hayvan sayısı
k <- nrow(cpped)       # Soyağacındaki hayvan sayısı

A <- as.matrix(nadiv::makeA(cpped[,1:3]))

# Sabit etkiler tasarım vektörleri
# Faktörler
xbay <- model.matrix(~-1+bay, data=cpsutver)
# Kovaryetler
xbyas <- as.numeric(cpsutver$byas)

# Sabit etkiler tasarım matrisi
X <- cbind(xbay, xbyas)
j <- ncol(X)               # Sabit etki sayısı

# Z matrisi
Z <- matrix(0, n, k)
for(i in 1:n){
   Z[i, cpped[i,1]] <- 1
}
# Fenotip vektörü
y <- as.vector(scale(cpsutver$sver305))

# Stan verisini oluştur
data_cpsutver <- list(N=n, J=j, K=k, X=X, Z=Z, Y=y, A=A)  

stanmodel1 <- stan(
  file = 'D:/lmmebook/stan/hayvanmodel.stan',  # Stan programı
  data = data_cpsutver,    # Stan verisi
  algorithm = "NUTS",      # Örnekleme algortiması
  chains = 4,              # Markov zincirleri sayısı
  iter = 40000,            # Zincir başına yineleme sayısı
  warmup = 20000,          # Zincir başına ısınma yinelemesi sayısı
  cores = 4                # İşlemci sayısı
)

summary(stanmodel1)

color_scheme_set("blue") 
rhats <- rhat(stanmodel1, pars=c("b", "h2", "sigma2_G", "sigma2_E"))
rhats
mcmc_rhat(rhats) + yaxis_text(hjust = 1)

h2_draws = rstan::extract(stanmodel1)$h2
# Sonsal ortalama (estimator)
mean(h2_draws)

# Sonsal aralıklar
quantile(h2_draws, probs=c(0.10, 0.90))

# Grafik
h2_draws_df <- data.frame(list(h2=h2_draws))
h2hist <- ggplot(h2_draws_df, aes(x=h2)) +
   geom_histogram(bins=20, color="gray")
h2hist



