library(rstan)
library(bayesplot)
library(loo)

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
  file="D:/lmmebook/stan/linearregression1.stan", 
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

# Doğruluk (Accuracy)
posteriors <- rstan::extract(stanmodel)
acc <- sum(round(posteriors$yhat)==round(y))/(4*n)
acc

## Stan’s optimization for estimation
# Max a posteriori (MAP)
# Max penalized likelihood (MLE).

stan_lm <- get_stancode(stanmodel1)
model = stan_model(model_code=stan_lm)
mle = optimizing(model, data=stan_data)
print(mle, digits=2)

# Örnek 4.4.15: Bir modeli yeniden kullanma
stanmodel2 <- stan(fit=stanmodel1, data=stan_data,
  iter=50, warmup=25, verbose=TRUE)

postmeans2 <- get_posterior_mean(stanmodel2)
head(postmeans2)
traceplot(stanmodel2)

# Parametreleri belli değerlerle başlatma
initpars <- function() {
  list(alpha = 50, 
  beta=array(c(2.1,2.5), dim=c(2,1)),
  sigma = 2)
}

stanmodel3 <- stan(fit=stanmodel1, data=stan_data,
  init=initpars,
  iter=50, warmup=25, verbose=TRUE)

postmeans3 <- get_posterior_mean(stanmodel3)
head(postmeans3)
traceplot(stanmodel3)

## Model Değerlendirme ve Karşılaştırma
##
library(loo)

# LOOCV
loo1 <- loo(stanmodel1, parameter_name="log_lik", save_psis = TRUE)
print(loo1)
plot(loo1, label_points = TRUE)

loo1a <- loo(stanmodel1, parameter_name="log_lik", is_method="sis", save_psis = TRUE)
print(loo1a)
plot(loo1a, label_points = TRUE)

loo2 <- loo(stanmodel2, parameter_name="log_lik", save_psis = TRUE)
print(loo2)
plot(loo2, label_points = TRUE)
loo3 <- loo(stanmodel3, parameter_name="log_lik", save_psis = TRUE)
print(loo3)
plot(loo3, label_points = TRUE)
loo_compare(loo1, loo2, loo3)

# WAIC
LL1  <- extract_log_lik(stanmodel1, parameter_name = "log_lik", merge_chains = TRUE)
waic1 <- waic(LL1)
print(waic1)
LL2  <- extract_log_lik(stanmodel2, parameter_name = "log_lik", merge_chains = TRUE)
waic2 <- waic(LL2)
print(waic2)
LL3  <- extract_log_lik(stanmodel3, parameter_name = "log_lik", merge_chains = TRUE)
waic3 <- waic(LL3)
print(waic3)
loo_compare(waic1, waic2, waic3)



## Denemeler

# S3 method for bayesQR
parameters::model_parameters(
  stanmodel,
  centrality = "median",
  dispersion = FALSE,
  ci = 0.95,
  ci_method = "eti",
  test = c("pd", "rope"),
  rope_range = "default",
  rope_ci = 0.95,
  #bf_prior = NULL,
  diagnostic = c("ESS", "Rhat"),
#  priors = TRUE,
  keep = NULL,
  drop = NULL,
  verbose = TRUE,
)

