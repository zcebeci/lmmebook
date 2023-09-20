
data {
  int<lower=0> N;       // Gözlem sayısı
  int<lower=0> K;       // Bağımsız değişken sayısı
  matrix[N, K] X;       // Tasarım matrisi
  real y[N];            // Gözlemler
}
parameters {
  real alpha;           // Kesme yüksekliği (intercept)
  vector[K] beta;       // Regresyon katsayıları
  real<lower=0> sigma;  // Hata standart sapması
}
model {
  vector[N] mu;
  mu <- alpha + X * beta;
  y ~ normal(mu, sigma); // Model (Hedef yoğunluk)
  // target += normal_lpdf(y | mu, sigma);
}
generated quantities{
 vector[N] log_lik;
 for (i in 1:N) log_lik[i] = normal_lpdf(y[i] | alpha + X[i, ] * beta, sigma);
}

