data {
  int<lower=0> N;       // Gözlem sayısı
  int<lower=0> K;       // Bağımsız değişken sayısı
  int<lower=0> S;       // Posterior sample size
  matrix[N, K] X;       // Tasarım matrisi (yeni veri)
  vector[S] alpha;      // Kesme yüksekliği (Bir önceki Stanmodel ile bulundu)
  matrix[S, K] beta;    // Regresyon katsayıları (Bir önceki Stanmodel ile bulundu)
}
parameters {
}
model {
}
generated quantities{
  matrix[S, N] y_test;
  for(i in 1:N) {
    for(j in 1:S) {
      y_test[j,i] = normal_rng(alpha[j] + X[i,] * beta[j,]);
    }  
  }
}

