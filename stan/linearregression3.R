# Bir Stan modeliyle yeni veri tahminleme
# Daha önceki bir Stan programından parametre tahminleri, tahmin programı için veri halinde iletilir.
# Yeni modelde parametre dağılımları tahmin edilmediği için linearregression2.stan kodunda
# parameters ve model öbeği boş bırakılır. Parametreler sabit olduğu için algorithm olarak "Fixed_param"
# kullanılır.
#
# Bu program linearregression1.R'dan sonra çalıştırılmalıdır.

library(rstan)
library(bayesplot)

# Test verisi (yeni veri)
set.seed(43)
n <-100  #Gözlem sayısı
x1 <- rnorm(n, 40, 10)  #1. bağımsız değişken 
x2 <- rnorm(n, 100, 20) #2. bağımsız değişken
Xnew <- cbind(x1, x2)

posteriors <- extract(stanmodel1)
alpha_post <- posteriors$alpha
beta_post <- posteriors$beta

stan_data <- list(
    N = nrow(Xnew),
    K = ncol(Xnew),
    S = length(alpha_post),
    X = Xnew,
    alpha=alpha_post,
    beta=beta_post
)
stantest <- rstan::stan(
  file="D:/lmmebook/stan/linearregression3.stan", 
  data=stan_data,
  algorithm = "Fixed_param"
)

print(stantest)
plot(stantest)
traceplot(stantest, pars=c("y_rep[1]", "y_rep[2]"))
postmeans <- get_posterior_mean(stantest)
head(postmeans)
tail(postmeans)


