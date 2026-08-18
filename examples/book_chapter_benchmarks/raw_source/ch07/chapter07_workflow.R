# =========================================================
# 第7章：联合物种分布模型的生物相互作用
# =========================================================

library(Hmsc)
# MCMC参数设置
nChains = 2
thin = 5
samples = 1000 
transient = 500*thin 
verbose = 500*thin



#模拟
n = 200
ns = 5
X = cbind(rep(1, n), rnorm(n), rnorm(n)) 
beta1 = rep(0, ns)
beta2 = c(2,2,-2,-2,0)
beta3 = c(1,-1,1,-1,0)
beta = cbind(beta1, beta2, beta3) 
L = X %*% t(beta)
Y = 1*((L + matrix(rnorm(n*ns), ncol = ns)) > 0)
colMeans(Y)
XData = data.frame(x1 = X[ ,2], x2 = X[ ,3]) 
studyDesign = data.frame(sample = as.factor(1:n)) 
rL = HmscRandomLevel(units = studyDesign$sample)



models = list() 
for (i in 1:3){
  XFormula = switch(i, ~1, ~x1, ~x1+x2)
  m = Hmsc(Y = Y, XData = XData, XFormula = XFormula, 
           studyDesign = studyDesign, 
           ranLevels = list(sample = rL),distr = "probit")
  models[[i]] = m
}
for (i in 1:3){
  models[[i]] = sampleMcmc(models[[i]], thin = thin,
                           samples = samples,transient = transient, 
                           nChains = nChains,verbose = verbose)
}
mpost = convertToCodaObject(models[[3]]) 
ess.beta = effectiveSize(mpost$Beta)
psrf.beta = gelman.diag(mpost$Beta, multivariate = FALSE)$psrf 
ess.omega = effectiveSize(mpost$Omega[[1]])
psrf.omega = gelman.diag(mpost$Omega[[1]], multivariate = FALSE)$psrf
library(corrplot)
for (i in 1:3){
  OmegaCor = computeAssociations(models[[i]]) 
  supportLevel = 0.95
  toPlot = ((OmegaCor[[1]]$support >  supportLevel)
            + (OmegaCor[[1]]$support < (1-supportLevel)) > 0)*OmegaCor[[1]]$mean
  corrplot(toPlot, method = "color", col= c("grey","white","black"))
}
partition = createPartition(m, nfolds = 2, column = "sample") 
partition.sp = c(1,2,3,4,5)
result = matrix(NA, nrow = 3, ncol = 3) 
for (i in 1:3){
  m = models[[i]]
  
  #Explanatory power
  preds = computePredictedValues(m)
  MF = evaluateModelFit(hM = m, predY = preds) 
  result[1,i] = mean(MF$TjurR2)
  
  #Predictive power based on cross-validation
  preds = computePredictedValues(m, partition = partition) 
  MF = evaluateModelFit(hM = m, predY = preds)
  result[2,i] = mean(MF$TjurR2)
  
  #Predictive power based on conditional cross-validation
  preds = computePredictedValues(m, partition = partition,
                                 partition.sp = partition.sp, mcmcStep = 100)
  MF = evaluateModelFit(hM = m, predY = preds) 
  result[3,i] = mean(MF$TjurR2)
}

#7.9 使用HMSC的真实案例研究：枯木栖息真菌的测序数据
data = read.csv("fungal data\\data.csv",stringsAsFactors=TRUE)
n = dim(data)[1]
head(data [, 1:6])
XData = data.frame(DC = as.factor(data$DC), 
                   readcount = data$readcount)
YData = data[,4:dim(data)[2]] 
sel.sp = colSums(YData > 0) >= 10 
YData = YData [, sel.sp]
P = colMeans(YData > 0)
A = colSums(YData)/sum(YData)
P = colMeans(YData>0)
A = colSums(YData)/sum(YData)

studyDesign = data.frame(sample = factor(data$LogID))
rL = HmscRandomLevel(units = studyDesign$sample) 
models = list()
for (i in 1:3){
  Y = as.matrix(YData)
  if (i==2) {Y = 1*(Y > 0)}
  if (i==3) {
    Y[Y==0] = NA
    Y = log(Y)
  }
  tmp = list() 
  for (j in 1:2){
    XFormula = switch(j, ~1 + log(readcount),~DC + log(readcount))
    m = Hmsc(Y = Y, XData = XData, XFormula = XFormula, 
             studyDesign = studyDesign, ranLevels = list(sample = rL), distr =  switch(i, "lognormal poisson", "probit", "normal"),YScale = TRUE)
    tmp[[j]] = m
  }
  models[[i]] = tmp
}
for (i in 1:3){
  for (j in 1:2){
    models[[i]][[j]] = sampleMcmc(models[[i]][[j]], thin = thin,
                                  samples = samples, transient = transient, nChains = nChains, 
                                  verbose = verbose, initPar = "ﬁxed effects")
  }
}
for (i in 1:3){
  mpost = convertToCodaObject(models[[i]][[2]]) 
  psrf.beta = gelman.diag(mpost$Beta,multivariate = FALSE)$psrf
  psrf.omega = gelman.diag(mpost$Omega[[1]], 
                           multivariate = FALSE)$psrf
}
m = models[[2]][[2]]
Gradient = constructGradient(m, focalVariable = "DC", 
                             non.focalVariables = list("readcount" = list(1)))
predY = predict(m, Gradient = Gradient, expected = TRUE)
for (i in 1:3){
  for (j in 1:2){
    OmegaCor = computeAssociations(models[[i]][[j]]) 
    supportLevel = 0.95
    toPlot = ((OmegaCor[[1]]$support > supportLevel)
              + (OmegaCor[[1]]$support < (1-supportLevel))> 0)*OmegaCor[[1]]$mean 
    corrplot(toPlot, method = "color", 
             col = c("grey","white","black"))
  }
}
for (j in 1:2){
  m = models[[2]][[j]]
  biPlot(m, etaPost = getPostEstimate(m, "Eta"), 
         lambdaPost = getPostEstimate(m, "Lambda"), colVar = 2)
}