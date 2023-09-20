# AlphaSimR paketinin kurulması/yüklenmesi
> suppressPackageStartupMessages(installed <- require(AlphaSimR))
> if (!installed) {
+   install.packages("AlphaSimR", repos="https://cloud.r-project.org")  
+   suppressPackageStartupMessages(require(AlphaSimR))
+ }

# Seleksiyon ve üreme simülasyonu
#
# Senaryo parametrelerini atama
> k <- 3 # Generasyon sayısı
> nSegs <- 100  # Kromozomlardak ayrılma bölgeleri sayısı
> mRate <- 1e-6 #Mutasyon oranı, NULL atanırsa uygulanmaz
> aRate <- NULL #Kayıp oranı, NULL atanırsa uygulanmaz
> fRate <- 0.8  #Seçilen dişilerin oranı
> mRate <- 0.1  #Seçilen erkeklerin oranı
> indTrait <- 1 # Özellik numarası
> nPro <- 1 # Çiftleşme başına yavru sayısı
> selCrit <- "bv"  #Seleksiyon kriteri: "gv","bv", "ebv", "pheno"
> popSize <- 100 # Populasyonda birey sayısı
> t1h2 <- 0.4   # Özellik 1 kalıtım derecesi
> t1Mean <- 500   # Özellik 1 ortalama
> t1Var <- 400   # Özellik 1 varyans

# Generasyonlar listesi oluşturma
> pop <- NULL
> pop <- vector(mode="list", length=k)

# Başlangıç/temel popülasyon
> print("Başlangıç popülasyonu oluşturuluyor.")
> set.seed(123)
> founderPop <- runMacs(
+   nInd = popSize,
+   nChr = 29,
+   segSites = nSegs,
+   inbred = FALSE,
+   species = "GENERIC",
+   split = 2,
+   ploidy = 2L,
+   manualCommand = NULL,
+   manualGenLen = NULL,
+   nThreads = NULL
+ )

# Simülasyon parametrelerini ayarla
> set.seed(123)
> SP <- SimParam$new(founderPop)
> SP$addTraitA(
+   nQtlPerChr = 5,
+   mean = t1Mean,
+   var = t1Var)
> SP$addSnpChip(nSnpPerChr=5)
> SP$traitNames <- (c("GCA"))
> SP$setSexes("yes_sys")
> SP$setVarE(h2 = t1h2)
# 1. popülasyonu oluştur
> pop[[1]] <- newPop(
+   founderPop,
+   simParam = SP)    

> nMale <- round(length(pop[[1]]@sex=="M")*mRate)
> nFemale <- round(length(pop[[1]]@sex=="F")*fRate)

# Seleksiyon ve çiftleştirme uygulaması
> print("Generasyonlar oluşturuluyor...")
> for (generation in 2:k) {
+    print(paste0("---", generation, ". generasyon oluşturuluyor."))
+    popSelM <- selectInd(
+      pop = pop[[generation-1]],
+      nInd = nMale,
+      trait = indTrait,
+      use = selCrit,
+      sex = "M",
+      selectTop = TRUE,
+      returnPop = TRUE,
+      candidates = NULL,
+      simParam = SP)
+   popSelF <- selectInd(
+      pop = pop[[generation-1]],
+      nInd = nFemale,
+      trait = indTrait,
+      use = selCrit,
+      sex = "F",
+      selectTop = TRUE,
+      returnPop = TRUE,
+      candidates = NULL,
+      simParam = SP)
+# Seçilenler popülasyonu oluşturma
+   popSel <- c(popSelM, popSelF)
+# Seçilenleri çiftleştirme
+   set.seed(123)
+   popPro <- randCross(
+      pop = popSel,
+      nCrosses = nFemale,
+      nProgeny = nPro,
+      simParam = SP)
+# Mutasyon uygulama
+   if(!is.null(mRate)){
+      popPro <- mutate(
+        pop = popPro, 
+        mutRate = mRate,
+        returnPos = FALSE, 
+        simParam = SP)
+   }
+# Doğum sonrası kayıplar
+   if(!is.null(aRate)){
+     popPro <- attrition(popPro, p = aRate)
+   }  
+   pop[[generation]] <- popPro
+   nMale <- round(length(pop[[generation]]@sex=="M")*mRate)
+   nFemale <- round(length(pop[[generation]]@sex=="F")*fRate)
+ }

# Populasyonları birleştir
> popList <- list(pop[[1]], pop[[2]]) 
> parentPop <- mergePops(popList)
> progenyPop <- pop[[k]]

# Parent pop genotip verisi
> parentGeno <- pullSegSiteGeno(parentPop, simParam=SP)
> progenyGeno <- pullSegSiteGeno(progenyPop, simParam=SP)

# Fenotip verisi 
> parentPheno <- pheno(parentPop)
> parentPheno <- data.frame(
+   id = parentPop@id,
+   sex = parentPop@sex,
+   Trait1 = parentPheno,
+   stringsAsFactors = FALSE)

> progenyPheno <- pheno(progenyPop)
> progenyPheno <- data.frame(
+   id = progenyPop@id,
+   sex = progenyPop@sex,
+   Trait1 = progenyPheno,
+   stringsAsFactors = FALSE)

# Soyağacı çekme
> parentPed <- getPed(parentPop)
> progenyPed <- getPed(progenyPop)

> wd <- getwd()
> setwd("D:/lmmebook/tinydata")
> write.table(parentGeno, file="CattleParentGeno.dat", 
+   sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)  
> write.table(parentPheno, file="CattleParentPheno.dat", 
+   sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)  
> write.table(parentPed, file="CattleParent.ped",
+   sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

> write.table(progenyGeno, file="CattleProgenyGeno.dat", 
+   sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)  
> write.table(progenyPheno, file="CattleProgenyPheno.dat", 
+   sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)  
> write.table(progenyPed, file="CattleProgeny.ped",
+   sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)


# PLINK dosyaları oluştur
> writePlink(
+   pop = parentPop,
+   baseName="CattleParentPlink",
+   traits = 1,
+   use = "pheno",
+   snpChip = 1,
+   useQtl = TRUE,
+   simParam = SP)

> writePlink(
+   pop = progenyPop,
+   baseName="CattleProgenyPlink",
+   traits = 1,
+   use = "pheno",
+   snpChip = 1,
+   useQtl = TRUE,
+   simParam = SP)

> setwd(wd)

