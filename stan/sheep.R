# Koyunlarda yapağı ağırlığı, lif çapı ve vücut ağırlığı
# Paketleri yükle
library(rstan)
library(nadiv)
library(bayesplot)
library(ggplot2)

# Fenotip dosyasını oku
sheep <- read.table("D:/lmmebook/datasets/sheep3.dat", header=TRUE)
sheep$animal <- as.factor(sheep$hayvan)
sheep$sire <- as.factor(sheep$sire)
sheep$dam <- as.factor(sheep$dam)
sheep$cins <- as.factor(sheep$cins)
sheep$dtip <- as.factor(sheep$dtip)
sheep$dyil <- as.factor(sheep$dyil)

# Soyağacı dosyasını oku
sheepped <- read.table("D:/lmmebook/datasets/sheepped3.dat", header=TRUE)
# Akrabalık matrisi
A <- as.matrix(nadiv::makeA(sheepped))

# X matrisi
xcins <- model.matrix(~-1+cins, data=sheep)
xdtip <- model.matrix(~-1+dtip, data=sheep)
xdyil <- model.matrix(~-1+dyil, data=sheep)
xvca <- sheep$vca
X <- cbind(xcins, xdtip, xdyil, xvca)

##################################################

n <- nrow(sheep)
k <- nrow(sheepped)
j <- ncol(X)

# Z matrisi
#Z <- model.matrix(~-1+as.factor(animal), data=sheep)
# Z matrisi
Z <- matrix(0, n, k)
for(i in 1:n)
   Z[i, k-n+i] <- 1

# Fenotipler
y1 <- as.vector(scale(sheep$typa))
y2 <- as.vector(scale(sheep$lcap))

# Stan'e gönderilecek veri
data_stan <- list(J=j, N=n, K=k, X=X, Z=Z, Y=y1, A=A) 

# Stan'i NUTS ile çalıştır

stanmodel1 <- stan(
  model_name <- "koyun_yapagi",
  file = 'D:/lmmebook/stan/hayvanmodel.stan', 
  data = data_stan, 
  chain = 4, 
  iter = 5000, 
  warmup = 2500,
  seed = 911)

summary(stanmodel1)
plot(stanmodel1, pars=c("h2", "sigma2_G", "sigma2_E"))

# Parametrelerin sonsal değerleri
h2_sonsal <- rstan::extract(stanmodel1)$h2
# Sonsal ortalama (estimator)
mean(h2_sonsal)

# Sonsal kantiller
quantile(h2_sonsal, probs=c(0.05, 0.5, 0.75, 0.95))

# Grafik
h2_sonsal_df <- data.frame(list(h2=h2_sonsal))
h2_hist <- ggplot(h2_sonsal_df, aes(x=h2)) +
   geom_histogram(bins=15, fill="dodgerblue", color="lightblue")
h2_hist

gtitle <- ggtitle("Sabit etkiler için sonsal dağılışlar
  (Medyan ve %90 interval")
mcmc_areas(stanmodel1, pars = c("h2"),
  prob = 0.90) + gtitle 

# Grafikler ve diyagnostik
color_scheme_set("blue") 
rhats <- rhat(stanmodel1, pars=c("b", "sigma2_G", "sigma2_E"))
rhats
mcmc_rhat(rhats) + yaxis_text(hjust = 1)


color_scheme_set("viridis")
mcmc_trace(
  as.array(stanmodel1),
  pars = c("h2",  "sigma_G","sigma_E"),
  size = 0.5,
  facet_args = list(ncol = 1),
  np = nuts_params(stanmodel1),
  np_style = trace_style_np(div_color = "red", div_size = 0.5)
)


color_scheme_set("green")
mcmc_hist(stanmodel1, pars = c("h2", "sigma_G","sigma_E"), binwidth=0.01)

color_scheme_set("blue")
mcmc_areas(stanmodel1, pars = c("h2", "sigma_G","sigma_E"))



