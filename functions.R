# Kod 3.1: A matrisi oluşturma fonksiyonu 1
buildA <- function(pedigree){
    if (nargs() > 1 ) stop("Soyağacı tablosu girilmeli")
    pedigree[is.na(pedigree)] <- 0
    baba <- pedigree[,2]
    ana <- pedigree[,3]
    n <- dim(pedigree)[1]
    m <- n + 1
    A <- matrix(0, ncol=m, nrow=m)
    baba <- (baba == 0)*(m) + baba
    ana <- (ana == 0)*m + ana
    for(i in 1:n){
       A[i,i] <- 1 + 0.5*A[baba[i], ana[i]]
       for (j in (i+1):n){    
         if (j > n) break
        A[j,i] <- A[i,j] <- 0.5*( A[i, baba[j]] + A[i, ana[j]] )
       }           
    }
    A <- A[-m, -m]
    return(A)  
}

# Kod 3.2: A matrisi oluşturma fonksiyonu 2
generateA <- function(pedigree){
   if (nargs() > 1 ) stop("Soyağacı tablosu girilmeli")
   pedigree[is.na(pedigree)] <- 0
   n <- nrow(pedigree)
   A <- matrix(0, n, n)
   F <- array(0, n)
   for(i in 1:n){
      A[i,i] <- 1
      if(pedigree[i, 2]!=0 && pedigree[i, 3]!=0){
         A[i, i] <- A[i, i] + 0.5*A[pedigree[i, 2], pedigree[i, 3]]
         A[1:i-1, i] <- 0.5*(A[1:i-1, pedigree[i,2]] + A[1:i-1, pedigree[i,3]])
      }else if(pedigree[i, 2]==0 && pedigree [i, 3]!=0){
         A[1:i-1, i] <- 0.5*(A[1:i-1, pedigree[i,3]])
      }else if(pedigree[i,2]!=0 && pedigree[i,3]==0){
        A[1:i-1, i] <- 0.5*(A[1:i-1, pedigree[i,2]])
      }
      A[i, 1:i-1] <- A[1:i-1, i]
      F[i] <- A[i,i]-1
   }
   return(list("A"=A, "F"=F))
} 

Kod 4.1: Farklı grafikleri karşılaştırma
> twoplots <- function(plot1, plot2, ncol = 2) {
    bayesplot::bayesplot_grid(
       plot1, plot2, 
       grid_args = list(ncol = ncol)
    )
}

Kod 4.2: Başlatma değerleri fonksiyonu
> initvals <- function(chain_id = 1) {
+  set.seed(chain_id)
+  list(
     Intercept = rnorm(n = 1, mean = 15, sd = 10),
     sigma = runif(n=1, min=0.4, max=0.8),
     sd = runif(n=1, min=0.3, max=0.6)
   )
}

# Kod 5.1:  Veri normalleştirme
rescale <- function(x, method="maxmin") {
   suppressPackageStartupMessages(installed <- require(clusterSim))
   if (!installed) {
     install.packages("clusterSim", repos="https://cloud.r-project.org")  
     suppressPackageStartupMessages(require(clusterSim))
   }
   if(method=="z"){
     xn <- clusterSim::data.Normalization(x, type="n1", normalization="column")
   }else if(method=="maxmin"){
     xn <- clusterSim::data.Normalization(x, type="n4", normalization="column")
     if(sum(is.na(xn))>0) stop("Veride Max=Min olan sütunlar var!")
   }else if(method=="max"){
     if(any(max(x)==0)) stop("Veride maksimumu 0 olan sütunlar var!")
     xn <- clusterSim::data.Normalization(x, type="n8", normalization="column")
   }else{
     stop("Yalnız max, maxmin ve z metotları kullanabilirsiniz!")
   }
   return(xn)
}


# Kod 5.2:  Normalleştirilmiş veriyi orijinale çevirme
back2original <- function(xn, method="maxmin") {
   unscale <- function(x)
     x*attr(xn,'normalized:scale')+attr(xn,'normalized:shift')
   if(method=="z"){
     xo <- t(apply(xn, 1, unscale))
   }else if(method=="maxmin"){
     xo <- t(apply(xn, 1, unscale))
   }else if(method=="max"){
     xo <- t(apply(xn, 1, unscale))
   }
   return(xo)
}

# Kod 5.3: Model başarım ölçütleri hesaplama
modelPerformance <- function(y, yhat, n, k){
   bias <- mean(y)-mean(yhat)
   errors <- y-yhat
   MAPE <- mean(abs(errors/y))
   REL <- cov(y,yhat)/sqrt(var(y))
   PCOR <- cor(y, yhat)
   R2 <- 1-(sum(errors^2)/sum((y-mean(y))^2))
   if(!missing(k))
     AdjR2 <- 1-((1-R2)*(n-1)/(n-k-1))
   else
     AdjR2 <- NA
   MAE <- mean(abs(errors))
   MAAPE <- mean(atan(abs(errors)/abs(y)))
   MSE <- mean(errors^2)
   RMSE <- sqrt(MSE)
   if(!missing(k)){
     AIC <- n * log(MSE) + 2 * k
     AICc <- AIC + ((2*k)*(k+1)/(n-k-1))
     BIC <- n * log(MSE) + k * log(n)
   }else{
     AIC <- AICc <- BIC <- NA
   }
   perf <- data.frame(Bias=bias, REL=REL, PCOR=PCOR, AdjR2=AdjR2,
     MAE=MAE, MAPE=MAPE, MAAPE=MAAPE, RMSE=RMSE, AIC=AIC, AICc=AICc, BIC=BIC)
     colnames(perf) <- c("Bias","REL", "PCOR", "AdjR2", "MAE", "MAPE", "MAAPE", 
       "RMSE", "AIC", "AICc", "BIC")
   return(perf)
}	

# Kod 5.4: Van Raden(2008) yöntemiyle G matrisi hesaplama 1
vanRadenG <- function(geno){ 
   pVar <- function(x){
     p <- mean(x)/2
     p <- 2*p*(1-p)
   }
   W <- scale(x=geno, center=TRUE, scale=FALSE)
   denom <- sum(apply(geno, 2, pVar))
   #denom <- sum(apply(X=W, 2, var)) #Alternatif hesaplama
   G <- tcrossprod(W)/denom
   return(G)
}

# Kod 5.5: Van Raden(2008) yöntemiyle G matrisi hesaplama 2
vanRadenG1 <- function(geno){ 
   n <- nrow(geno)
   m <- ncol(geno)
   p <- colMeans(geno)/2 #MAF
   denom <- 0
   for(j in 1:m)
     denom <- denom + 2*p[j]*(1-p[j])
   W <- matrix(0, nrow = n, ncol = m)
   for(i in 1:n){
     for(j in 1:m){
       W[i,j] <- geno[i,j] - 2*p[j]
     }
   }
   G <- tcrossprod(W)/denom
   return(G)
}

# Kod 5.6: Van Raden(2008) 2. yöntemiyle G matrisi hesaplama
vanRadenG2 <- function(geno){
   Z <- scale(x=as.matrix(geno), center=TRUE, scale=TRUE)
   if(length(which(is.na(Z)))>0){
     Z[which(is.na(Z))] <- 0
     cat("There are NAs in Z! They were assumed 0. Use another method for G")
   }
   G <- tcrossprod(Z)/ncol(Z)
   return(G)
}

# Kod 5.7: Forni(2011) yöntemiyle G matrisi
forniG <- function(geno, option=1){
   m <- ncol(geno)
   if(option==1){ 
     p <- 0.5
     W <- scale(x=geno, center=rep(2*p, m), scale=FALSE)
     Dnom <- 2*p*(1-p)*m 
     G <- tcrossprod(W)/Dnom
   }else if(option==2){ 
     p <- apply(X=geno, 2, FUN=function(x){p=mean(x)/2})
     p <- ifelse(p>0.5, (1-p), p)
     p <- mean(p)
     W <- scale(x=geno, center=rep(2*p, m), scale=FALSE)
     Dnom <- 2*p*(1-p)*m 
     G <- tcrossprod(W)/Dnom
   }else if(option==3){
     W <- scale(x=geno, center=TRUE, scale=FALSE)
     WW <- tcrossprod(W)
     Dnom <- sum(diag(WW))/m
     G <- G/Dnom
   }
   return(G)
}

# Kod 5.8: GBLUP fonksiyonu
gblup <- function(y, X, Z, G, varG, varE) {
   XX <- crossprod(X)		
   XZ <- crossprod(X, Z)
   ZX <- crossprod(Z, X)
   ZZ <- crossprod(Z)
   Xy <- crossprod(X, y)
   Zy <- crossprod(Z, y)
   lambda <- varE/varG
   Ginv <- solve(G)
   ZZGinv <- ZZ + (Ginv * lambda)
   LHS1 <- cbind(XX, XZ)
   LHS2 <- cbind(ZX, ZZGinv)
   LHS <- rbind(LHS1, LHS2)
   RHS <- rbind(Xy, Zy)
   sol <- solve(LHS, RHS)
   return(sol)
}

# Kod 5.9: Çapraz Doğrulama
# Bağımlılık – Kod 5.3
#
cvGBLUP <-function(y, X, W, tRate=0.7, niter=5){
   n <- nrow(W)
   nTrain <- round(n*tRate)
   metrics <- as.data.frame(matrix(NA, nrow=niter+1, ncol=10))
   colnames(metrics) <- c("Bias", "REL", "PCOR", "R2",
     "MAE", "MAAPE", "RMSE", "AIC", "AICc", "BIC")
   rownames(metrics) <- c(paste0("Iter", 1:niter), "Ort.")
   for (i in 1:niter){
     trainIdx <- sample(1:n, size=nTrain, replace=FALSE)
     trainY <- y[trainIdx]
     trainW <- W[trainIdx,]
     trainX <- X[trainIdx,]
     testY <- y[-trainIdx]
     testW <- W[-trainIdx,]
     testX <- X[-trainIdx,]
     trainG <- vanRadenG(trainW)
     trainZ <- diag(nrow(trainW))
     testZ <- diag(nrow(testW))
     Z <- diag(nrow(W))
     G <- vanRadenG(W)
     modelGBLUP <- mixed.solve(
       method="REML",
       y = trainY, 
       X = trainX,
       Z = trainZ,
       K = trainG, 
       SE = TRUE, 
       return.Hinv = FALSE)
     trainYhat <-  
       tcrossprod(modelGBLUP$u,trainZ) + tcrossprod(modelGBLUP$beta,trainX)
     k <- length(modelGBLUP$u) + length(modelGBLUP$beta)
     metrics[i,] <- modelPerformance(trainY, as.vector(trainYhat), n,k)
   }
   metrics[niter+1,] <- apply(as.matrix(metrics), 2, mean, na.rm=T)
   return(metrics)
}

