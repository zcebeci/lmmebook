data {
  int<lower=1>    J; // Sabit etkilerin sayısı
  int<lower=1>    K; // Soyağacındaki hayvanların sayısı
  int<lower=1>    N; // Fenotip (gözlem) sayısı
  matrix[N,J]     X; // Sabit etkiler için tasarım matrisi
  matrix[N,K]     Z; // Rastlantısal etkiler için tasarım matrisi
  matrix[K,K]     A; // Akrabalık matrisi
  vector[N]       Y; // Fenotip (Bağımlı değişken) 
}
transformed data{
  matrix[K,K] LA;
  LA = cholesky_decompose(A);
}
parameters {
  vector[K]  u; // Damızlık değerleri
  vector[J] b; // Sabit etkiler
  real<lower=0> sigma_G; // Genetik standart sapma
  real<lower=0> sigma_E; // Hata standart sapması
}
model {
    vector[N] mu;
    vector[K] a;
    u ~ normal(0, 1);
    a = sigma_G * (LA * u);
    mu = X * b + Z * a;
    Y ~ normal(mu, sigma_E);
    to_vector(b) ~ normal(0, 1);
    sigma_G ~ student_t(4, 0, 1);
    sigma_E ~ student_t(4, 0, 1);
}
generated quantities{
  real sigma2_G;
  real sigma2_E;
  real h2;
  sigma2_G = sigma_G * sigma_G; // Genetik varyans
  sigma2_E = sigma_E * sigma_E; // Hata varyansı
  h2 = sigma2_G / (sigma2_G + sigma2_E); // Kalıtım derecesi
  vector[N] log_lik;
  for (i in 1:N) log_lik[i] = normal_lpdf(Y[i] | X[i, ] * b + Z[i,] * u, sigma_E);
}
