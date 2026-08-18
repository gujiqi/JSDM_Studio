# Set the base directory using your favorite method
# setwd("...")

##################################################################################################
# INPUT AND OUTPUT OF THIS SCRIPT (BEGINNING)
##################################################################################################
#	INPUT. Fitted models

#	OUTPUT. MCMC convergence statistics for selected model parameters,
# illustrated (for all RUNs performed thus far in S3) in the file "results/MCMC_convergence.pdf",
# and the text file "results/MCMC_convergence.txt".
##################################################################################################
# INPUT AND OUTPUT OF THIS SCRIPT (END)
##################################################################################################


##################################################################################################
# MAKE THE SCRIPT REPRODUCIBLE (BEGINNING)
##################################################################################################
set.seed(1)
##################################################################################################
## MAKE THE SCRIPT REPRODUCIBLE (END)
##################################################################################################
localDir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
##################################################################################################
# MAKE THE SCRIPT REPRODUCIBLE (BEGINNING)
##################################################################################################
set.seed(1)
##################################################################################################
## MAKE THE SCRIPT REPRODUCIBLE (END)
##################################################################################################

load("models/models_thin_10_samples_1000_chains_4.Rdata")
names(models)

##################################################################################################
# SETTING COMMONLY ADJUSTED PARAMETERS TO NULL WHICH CORRESPONDS TO DEFAULT CHOICE (BEGINNING)
##################################################################################################
showBeta = NULL #Default: showBeta = TRUE, convergence shown for beta-parameters
showGamma = NULL #Default: showGamma = FALSE, convergence not shown for gamma-parameters
showOmega = NULL #Default: showOmega = FALSE, convergence not shown for Omega-parameters
maxOmega = NULL #Default: convergence of Omega shown for 50 randomly selected species pairs
showRho = NULL #Default: showRho = FALSE, convergence not shown for rho-parameters
showAlpha = NULL #Default: showAlpha = FALSE, convergence not shown for alpha-parameters
##################################################################################################
# SETTING COMMONLY ADJUSTED PARAMETERS TO NULL WHICH CORRESPONDS TO DEFAULT CHOICE (END)
##################################################################################################

##################################################################################################
# CHANGE DEFAULT OPTIONS BY REMOVING COMMENT AND SETTING VALUE (BEGINNING)
# NOTE THAT THIS IS THE ONLY SECTION OF THE SCRIPT THAT YOU TYPICALLY NEED TO MODIFY
##################################################################################################
showBeta = TRUE
showGamma = TRUE
showOmega = TRUE
maxOmega = 100
showRho = TRUE
showAlpha = TRUE
##################################################################################################
# CHANGE DEFAULT OPTIONS BY REMOVING COMMENT AND SETTING VALUE (END)
# NOTE THAT THIS IS THE ONLY SECTION OF THE SCRIPT THAT YOU TYPICALLY NEED TO MODIFY
##################################################################################################

##################################################################################################
# SET DIRECTORIES (BEGINNING)
##################################################################################################
localDir = "."
modelDir = file.path(localDir, "models")
resultDir = file.path(localDir, "results")
if (!dir.exists(resultDir)) dir.create(resultDir)
##################################################################################################
# SET DIRECTORIES (END)
##################################################################################################

if(is.null(showBeta)) showBeta = TRUE
if(is.null(showGamma)) showGamma = FALSE
if(is.null(showOmega)) showOmega = FALSE
if(is.null(maxOmega)) maxOmega = 50
if(is.null(showRho)) showRho = FALSE
if(is.null(showAlpha)) showAlpha = FALSE

library(Hmsc)
library(colorspace)
library(vioplot)

samples_list = c(5,250,250,250,250,250)
thin_list = c(1,1,10,100,1000,10000)
nst = length(thin_list)
nChains = 4

text.file = file.path(resultDir,"/MCMC_convergence.txt")
cat("MCMC Convergennce statistics\n\n",file=text.file,sep="")

ma.beta = NULL
na.beta = NULL
ma.gamma = NULL
na.gamma = NULL
ma.omega= NULL
na.omega = NULL
ma.alpha = NULL
na.alpha = NULL  
ma.rho = NULL
na.rho = NULL
Lst = 1
while(Lst <= nst){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  
  
  filename = file.path(modelDir,paste("models_thin_", as.character(thin),
                                      "_samples_", as.character(samples),
                                      "_chains_",as.character(nChains),
                                      ".Rdata",sep = ""))
  if(file.exists(filename)){
    load(filename)
    cat(c("\n",filename,"\n\n"),file=text.file,sep="",append=TRUE)
    nm = length(models)
    for(j in 1:nm){
      mpost = convertToCodaObject(models[[j]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
      nr = models[[j]]$nr
      cat(c("\n",names(models)[j],"\n\n"),file=text.file,sep="",append=TRUE)
      if(showBeta){
        psrf = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
        tmp = summary(psrf)
        cat("\nbeta\n\n",file=text.file,sep="",append=TRUE)
        cat(tmp[,1],file=text.file,sep="\n",append=TRUE)
        if(is.null(ma.beta)){
          ma.beta = psrf[,1]
          na.beta = paste0(as.character(thin),",",as.character(samples))
        } else {
          ma.beta = cbind(ma.beta,psrf[,1])
          if(j==1){
            na.beta = c(na.beta,paste0(as.character(thin),",",as.character(samples)))
          } else {
            na.beta = c(na.beta,"")
          }
        }
      }
      if(showGamma){
        psrf = gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf
        tmp = summary(psrf)
        cat("\ngamma\n\n",file=text.file,sep="",append=TRUE)
        cat(tmp[,1],file=text.file,sep="\n",append=TRUE)
        if(is.null(ma.gamma)){
          ma.gamma = psrf[,1]
          na.gamma = paste0(as.character(thin),",",as.character(samples))
        } else {
          ma.gamma = cbind(ma.gamma,psrf[,1])
          if(j==1){
            na.gamma = c(na.gamma,paste0(as.character(thin),",",as.character(samples)))
          } else {
            na.gamma = c(na.gamma,"")
          }
        }
      }
      if(showRho & !is.null(mpost$Rho)){
        psrf = gelman.diag(mpost$Rho,multivariate=FALSE)$psrf
        cat("\nrho\n\n",file=text.file,sep="",append=TRUE)
        cat(psrf[1],file=text.file,sep="\n",append=TRUE)
      }
      if(showOmega & nr>0){
        cat("\nomega\n\n",file=text.file,sep="",append=TRUE)
        for(k in 1:nr){
          cat(c("\n",names(models[[j]]$ranLevels)[k],"\n\n"),file=text.file,sep="",append=TRUE)
          tmp = mpost$Omega[[k]]
          z = dim(tmp[[1]])[2]
          if(z > maxOmega){
            sel = sample(1:z, size = maxOmega)
            for(i in 1:length(tmp)){
              tmp[[i]] = tmp[[i]][,sel]
            }
          }
          psrf = gelman.diag(tmp, multivariate = FALSE)$psrf
          tmp = summary(psrf)
          cat(tmp[,1],file=text.file,sep="\n",append=TRUE)
          if(is.null(ma.omega)){
            ma.omega = psrf[,1]
            na.omega = paste0(as.character(thin),",",as.character(samples))
          } else {
            ma.omega = cbind(ma.omega,psrf[,1])
            if(j==1){
              na.omega = c(na.omega,paste0(as.character(thin),",",as.character(samples)))
            } else {
              na.omega = c(na.omega,"")
            }
          }
        }
      }
      if(showAlpha & nr>0){
        for(k in 1:nr){
          if(models[[j]]$ranLevels[[k]]$sDim>0){
            cat("\nalpha\n\n",file=text.file,sep="\n",append=TRUE)
            cat(c("\n",names(models[[j]]$ranLevels)[k],"\n\n"),file=text.file,sep="",append=TRUE)
            psrf = gelman.diag(mpost$Alpha[[k]],multivariate = FALSE)$psrf
            cat(psrf[,1],file=text.file,sep="\n",append=TRUE)            
          }
        }
      }
    }
  }
  Lst = Lst + 1
}

pdf(file= file.path(resultDir,"/MCMC_convergence.pdf"))
if(showBeta){
  par(mfrow=c(2,1))
  vioplot(ma.beta,col=rainbow_hcl(nm),names=na.beta,ylim=c(0,max(ma.beta)),main="psrf(beta)")
  legend("topright",legend = names(models), fill=rainbow_hcl(nm))
  vioplot(ma.beta,col=rainbow_hcl(nm),names=na.beta,ylim=c(0.9,1.1),main="psrf(beta)")
}
if(showGamma){
  par(mfrow=c(2,1))
  vioplot(ma.gamma,col=rainbow_hcl(nm),names=na.gamma,ylim=c(0,max(ma.gamma)),main="psrf(gamma)")
  legend("topright",legend = names(models), fill=rainbow_hcl(nm))
  vioplot(ma.gamma,col=rainbow_hcl(nm),names=na.gamma,ylim=c(0.9,1.1),main="psrf(gamma)")
}
if(showOmega & !is.null(ma.omega)){
  par(mfrow=c(2,1))
  vioplot(ma.omega,col=rainbow_hcl(nm),names=na.omega,ylim=c(0,max(ma.omega)),main="psrf(omega)")
  legend("topright",legend = names(models), fill=rainbow_hcl(nm))
  vioplot(ma.omega,col=rainbow_hcl(nm),names=na.omega,ylim=c(0.9,1.1),main="psrf(omega)")
}
dev.off()












##########################################################################################################




################################################################################
# Hmsc-HPC MCMC convergence diagnostics (PSRF / Gelman-Rubin)
# Adapted for fitted models imported from HPC
# Author: YOU
################################################################################

##################################################################################################
# 1. MAKE THE SCRIPT REPRODUCIBLE
##################################################################################################
set.seed(1)

##################################################################################################
# 2. USER-DEFINED SETTINGS  ⭐⭐⭐（你只需要改这里）
##################################################################################################

# 你真实使用的 MCMC 参数
samples_list <- c(1000)
thin_list    <- c(10)
nChains      <- 4

# 只画 Beta（强烈推荐）
showBeta  <- TRUE
showGamma <- FALSE
showOmega <- FALSE
showRho   <- FALSE
showAlpha <- FALSE
maxOmega  <- 50

##################################################################################################
# 3. SET DIRECTORIES
##################################################################################################
localDir  <- "."
modelDir  <- file.path(localDir, "models")
resultDir <- file.path(localDir, "results")

if (!dir.exists(modelDir)) {
  stop("❌ models 目录不存在，请检查路径")
}
if (!dir.exists(resultDir)) dir.create(resultDir)

##################################################################################################
# 4. LOAD PACKAGES
##################################################################################################
library(Hmsc)
library(colorspace)
library(vioplot)

##################################################################################################
# 5. INITIALIZE OUTPUT
##################################################################################################
text.file <- file.path(resultDir, "MCMC_convergence.txt")
cat("MCMC convergence diagnostics (PSRF)\n\n", file = text.file)

ma.beta <- NULL
na.beta <- NULL

##################################################################################################
# 6. LOOP OVER MCMC SETTINGS
##################################################################################################
nst <- length(thin_list)

for (Lst in seq_len(nst)) {
  
  thin    <- thin_list[Lst]
  samples <- samples_list[Lst]
  
  filename <- file.path(
    modelDir,
    paste0(
      "models_thin_", thin,
      "_samples_", samples,
      "_chains_", nChains,
      ".Rdata"
    )
  )
  
  if (!file.exists(filename)) {
    warning("⚠️ 找不到文件：", filename)
    next
  }
  
  load(filename)   # 必须加载出对象：models
  nm <- length(models)
  
  cat("\n====================================================\n",
      filename, "\n", file = text.file, append = TRUE)
  
  for (j in seq_len(nm)) {
    
    cat("\nModel: ", names(models)[j], "\n",
        file = text.file, append = TRUE)
    
    mpost <- convertToCodaObject(
      models[[j]],
      spNamesNumbers  = c(TRUE, FALSE),
      covNamesNumbers = c(TRUE, FALSE)
    )
    
    if (showBeta) {
      psrf <- gelman.diag(mpost$Beta, multivariate = FALSE)$psrf[, 1]
      
      cat("\nBeta PSRF summary:\n",
          file = text.file, append = TRUE)
      cat(summary(psrf), file = text.file, append = TRUE)
      
      if (is.null(ma.beta)) {
        ma.beta <- psrf
        na.beta <- paste0(thin, ",", samples)
      } else {
        ma.beta <- cbind(ma.beta, psrf)
        na.beta <- c(na.beta, paste0(thin, ",", samples))
      }
    }
  }
}

##################################################################################################
# 7. SAFETY CHECK
##################################################################################################
if (is.null(ma.beta)) {
  stop("❌ 没有成功读入任何模型，请检查 models_thin_* 文件名是否匹配")
}

##################################################################################################
# 8. PLOT PSRF (BETA ONLY)
##################################################################################################
pdf(file = file.path(resultDir, "MCMC_convergence.pdf"), width = 8, height = 6)

par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

# 全范围
vioplot(
  ma.beta,
  col   = rainbow_hcl(nm),
  names = na.beta,
  ylim  = c(0, max(ma.beta)),
  main  = "PSRF (Beta) - full range",
  ylab  = "PSRF"
)
legend("topright", legend = names(models), fill = rainbow_hcl(nm))

# 放大到 0.9–1.1（论文最关键）
vioplot(
  ma.beta,
  col   = rainbow_hcl(nm),
  names = na.beta,
  ylim  = c(0.9, 1.1),
  main  = "PSRF (Beta) - zoomed",
  ylab  = "PSRF"
)
abline(h = 1.05, col = "red", lty = 2)

dev.off()

cat("\n✅ MCMC convergence diagnostics finished.\n")







