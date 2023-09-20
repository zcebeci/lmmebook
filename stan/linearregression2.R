library(rstan)
library(bayesplot)

# Benzetimle veri oluşturma
set.seed(123)
n <-100  #Gözlem sayısı
x1 <- rnorm(n, 40, 10)  #1. bağımsız değişken 
x2 <- rnorm(n, 100, 20) #2. bağımsız değişken
e <- rnorm(n, 5, 2)  # Rastlantısal hata
a <- 50  # Kesme yüksekliği (intercept)
b1 <- 3  # x1 için regresyon katsayısı
b2 <- 2  # x2 için regresyon katsayısı
y <- a + b1*x1 + b2*x2 + e  # Bağımlı değişken
df <- data.frame(y, x1, x2)  # Veri çerçevesi
pairs(df)  # İlişki grafikleri

X <- df[, c("x1","x2")]
y <- df[, "y"]

# Verinin hazırlanması
stan_data <- list(
    N = nrow(X),
    K = ncol(X),
    X = X,
    y = y
)

# Stan'in çalıştırılması
stanmodel1 <- rstan::stan(
  file="D:/lmmebook/stan/linearregression2.stan", 
  data=stan_data
)

# Sonuçların incelenmesi
print(stanmodel1)
plot(stanmodel1)
traceplot(stanmodel1)
plot(stanmodel1, pars = c("alpha", "beta"))
traceplot(stanmodel1, pars = c("alpha", "beta"))
postmeans <- get_posterior_mean(stanmodel1)
head(postmeans)
tail(postmeans)

color_scheme_set("green")
y_rep <- as.matrix(stanmodel1, pars = "y_rep")
dim(y_rep)
ppc_dens_overlay(y, y_rep[1:200, ])

color_scheme_set("brightblue")
bayesplot::ppc_stat(y = y, yrep = y_rep, stat = "mean")
bayesplot::ppc_stat(y = y, yrep = y_rep, stat = "var")

color_scheme_set("purple")
bayesplot::ppc_scatter_avg(y = y, yrep = y_rep)
