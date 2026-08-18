# =========================================================
# 第6章：联合物种分布模型：物种生态位变异
# =========================================================

# 加载所需R包
library(ape)      # 用于系统发育树的生成和处理
library(MASS)     # 提供多元正态分布函数 mvrnorm
library(Hmsc)     # Hierarchical Modelling of Species Communities 包

# 设置物种数量
ns = 100

# 生成一个共祖模型的系统发育树（100个物种）
phy = rcoal(n = ns, tip.label = 
              sprintf('sp_%.3d',1:ns), br = "coalescent")

# 根据Brownian motion模型计算系统发育协方差矩阵
C = vcv(phy, model = "Brownian", corr = TRUE)

# 构造两种不同的物种性状结构
# A: 性状与系统发育无关（随机）
Tr.A = cbind(rep(1,ns), rnorm(ns))

# B: 性状与系统发育相关（通过协方差矩阵模拟）
Tr.B = cbind(rep(1,ns), mvrnorm(n = 1,
                                mu = rep(0, ns), Sigma = C))

# 物种性状对环境响应的参数（gamma）
gamma = cbind(c(-2,2), c(-1,1))

# 计算环境响应的均值
mu.A = gamma %*% t(Tr.A)
mu.B = gamma %*% t(Tr.B)

# 环境响应的协方差矩阵
V2 = diag(2)

# 模拟物种对环境的回归系数
# A: 系数具有系统发育结构
beta.A = matrix(mvrnorm(n = 1, mu = as.vector(mu.A), 
                        Sigma = kronecker(C, V2)), ncol = ns)

# B: 系数不具有系统发育结构
beta.B = matrix(mvrnorm(n = 1, mu = as.vector(mu.B), 
                        Sigma = kronecker(diag(ns), V2)), ncol = ns)

# 设置样点数量
n = 50

# 构造环境变量矩阵
X = cbind(rep(1, n), rnorm(n))

# 计算潜在生态位值（线性预测）
L.A = X %*% beta.A
L.B = X %*% beta.B

# 根据潜在变量生成物种存在/缺失数据（probit模型）
Y.A = 1*((L.A + matrix(rnorm(n*ns), ncol = ns)) > 0)
Y.B = 1*((L.B + matrix(rnorm(n*ns), ncol = ns)) > 0)

# 计算每个样点的物种丰富度
S.A = rowSums(Y.A)

# 计算每个物种的出现概率
P.A = colMeans(Y.A)

S.B = rowSums(Y.B)
P.B = colMeans(Y.B)

# 选择模拟的群落类型
community = "A"

# 根据群落类型选择对应数据
Y = switch(community, "A" = Y.A, "B" = Y.B)

# 设置物种名称
colnames(Y) = phy$tip.label

# 选择对应的物种性状
Tr = switch(community, "A" = Tr.A, "B" = Tr.B)

# 构建物种性状数据框
TrData = data.frame(trait = Tr [,2])
rownames(TrData) = phy$tip.label

# 构建环境变量数据框
XData = data.frame(x = X [,2])

# 构建HMSC模型
m= Hmsc(Y = Y, XData= XData, XFormula = ~x, TrData = TrData,
        TrFormula = ~trait, phyloTree = phy, distr = "probit")

# MCMC参数设置
nChains = 2
thin = 5
samples = 1000 
transient = 500*thin 
verbose = 500*thin

# 运行MCMC采样
m = sampleMcmc(m, thin = thin, samples = samples, 
               transient = transient, nChains = nChains, verbose = verbose)

# 转换为coda对象以便进行MCMC诊断
mpost =  convertToCodaObject(m)

# 计算有效样本量
effectiveSize(mpost$Rho)

# Gelman-Rubin收敛诊断
gelman.diag(mpost$Rho,multivariate=FALSE, autoburnin=FALSE) $psrf

# 计算模型预测值
preds = computePredictedValues(m)

# 评估模型拟合优度
MF = evaluateModelFit(hM = m, predY = preds)

# 创建交叉验证分区
partition = createPartition(m, nfolds = 2)

# 计算交叉验证预测
preds = computePredictedValues(m, partition = partition)

# 评估交叉验证模型性能
MFCV = evaluateModelFit(hM = m, predY = preds)

# 提取Beta参数的后验估计
postBeta = getPostEstimate(m, parName = "Beta")

# 绘制物种环境响应系数（带系统发育树）
plotBeta(m, post = postBeta, param = "Sign", plotTree = TRUE,
         supportLevel = 0.95, split = 0.4, spNamesNumbers = c(F,F))

# 提取Gamma参数（性状与环境关系）
postGamma = getPostEstimate(m, parName = "Gamma")

# 绘制Gamma结果
plotGamma(m, post = postGamma, param = "Sign",
          supportLevel = 0.95)

# 查看系统发育信号rho的后验分布
summary(mpost$Rho) $quantiles

# 构建不包含性状和系统发育信息的模型（Null trait phylogeny model）
m = Hmsc(Y = Y, XData = XData, XFormula = ~x, distr = "probit")

# 重新运行MCMC
m = sampleMcmc(m, thin = thin, samples = samples, transient = transient,
               nChains = nChains, verbose = verbose)

# 计算模型预测
preds = computePredictedValues(m)

# 计算模型拟合
MF.NTP = evaluateModelFit(hM = m, predY = preds)

# 交叉验证预测
preds = computePredictedValues(m, partition = partition) 

# 交叉验证模型性能
MFCV.NTP = evaluateModelFit(hM = m, predY = preds)

# 比较包含性状+系统发育 与 不包含模型 的预测能力差异
Delta.TjurR2 = MFCV$TjurR2-MFCV.NTP$TjurR2

# 提取真实斜率参数
beta.slope.true = beta.A[2,]

# 提取估计斜率
beta.slope.est = postBeta$mean [2,]

# 提取无系统发育模型的Beta估计
postBeta.NTP = getPostEstimate(m, parName = "Beta")

# 提取其斜率
beta.slope.est.NTP = postBeta.NTP$mean[2,]

# 切换为群落 B
community = "B"

# 重新选择群落数据
Y = switch(community, "A" = Y.A, "B" = Y.B)

# 设置物种名称
colnames(Y) = phy$tip.label

# 选择对应的物种性状
Tr = switch(community, "A" = Tr.A, "B" = Tr.B)

# 构建物种性状数据
TrData = data.frame(trait = Tr[,2])
rownames(TrData) = phy$tip.label

# 环境变量数据（不需要改变）
XData = data.frame(x = X[,2])

# 重新构建 HMSC 模型
m= Hmsc(Y = Y, XData= XData, XFormula = ~x, TrData = TrData,
        TrFormula = ~trait, phyloTree = phy, distr = "probit")

# 再次运行MCMC
m = sampleMcmc(m, thin = thin, samples = samples,
               transient = transient,
               nChains = nChains, verbose = verbose)

# 转换为coda对象
mpost = convertToCodaObject(m)

# 提取Beta参数
postBeta = getPostEstimate(m, parName = "Beta")

# 绘制Beta
plotBeta(m, post = postBeta, param = "Sign",
         plotTree = TRUE, supportLevel = 0.95,
         split = 0.4, spNamesNumbers = c(F,F))

# 提取Gamma
postGamma = getPostEstimate(m, parName = "Gamma")

# 绘制Gamma
plotGamma(m, post = postGamma, param = "Sign",
          supportLevel = 0.95)

# 查看系统发育信号
summary(mpost$Rho)$quantiles

# 读取真实植物数据
data = read.csv("plant data\\whittaker revisit data.csv",stringsAsFactors = TRUE)

# 查看数据前几行
head(data)

# 将site转换为因子
data$site = factor(data$site)

# 提取所有site
sites = levels(data$site)

# 提取所有物种
species = levels(data$species)

# site数量
n = length(sites)

# 物种数量
ns = length(species)

# 初始化物种矩阵
Y = matrix(NA, nrow = n, ncol = ns)

# 初始化环境变量
env = rep(NA, n)

# 初始化性状
trait = rep(NA, ns)

# 构建site × species矩阵
for (i in 1:n){
  for (j in 1:ns){
    row = data$site==sites[i] & data$species==
      species[j]
    Y[i,j] = data[row,]$value
    env[i] = data[row,]$env
    trait[j] = data[row,]$trait
  }
}

# 设置物种名称
colnames(Y) = species

# 构建环境变量数据
XData = data.frame(TMG = env)

# 构建性状数据
TrData = data.frame(CN = trait, row.names = species)

# 计算物种出现概率
P = colMeans(Y > 0)

# 计算平均丰度
A = colSums(Y)/colSums(Y > 0)

# 读取分类学信息
taxonomy = read.csv("plant data\\taxonomy.csv") 

# 构建系统发育树
library(ape)
# 转换成分类变量。
taxonomy$family <- factor(taxonomy$family)
taxonomy$genus <- factor(taxonomy$genus)
taxonomy$species <- factor(taxonomy$species)
plant.tree <- as.phylo(~family/genus/species, data=taxonomy, collapse = FALSE)
plant.tree$edge.length <- rep(1,length(plant.tree$edge))

# 设置模型公式
XFormula = ~TMG
TrFormula = ~CN

# 创建模型列表
models = list()

# 存在/缺失模型
models[[1]] = Hmsc(Y=1*(Y > 0), XData = XData,
                   XFormula = XFormula,TrData = TrData,
                   TrFormula = TrFormula, phyloTree = plant.tree,
                   distr = "probit")

# 丰度模型
models[[2]] = Hmsc(Y = Y, XData = XData,
                   XFormula = XFormula,TrData = TrData,
                   TrFormula = TrFormula, phyloTree = plant.tree,
                   distr = "lognormal poisson")

# 运行两个模型的MCMC
for (i in 1:2){
  models[[i]] = sampleMcmc(models[[i]], thin = thin, 
                           samples = samples, transient = transient,
                           nChains = nChains, verbose = verbose)
}

# 沿环境梯度进行预测
for (i in 1:2){
  m = models [[i]]
  Gradient = constructGradient(m, focalVariable = "TMG")
  predY = predict(m, Gradient = Gradient, expected = TRUE)
  q = c(0.25,0.5,0.75)
  
  # 绘制物种丰富度沿梯度变化
  plotGradient(m, Gradient, pred = predY, measure = "S", 
               showData = TRUE, q = q)
  
  # 绘制单个物种响应
  plotGradient(m, Gradient, pred = predY, measure
               = "T", index = 2,showData = TRUE, q = q)
}

# 方差分解
VP = computeVariancePartitioning(models [[1]],
                                 group = c(1,1),   groupnames = "TMG")

# 查看解释度
VP$R2T

VP = computeVariancePartitioning(models[[2]],
                                 group = c(1,1),groupnames = "TMG")

VP$R2T

# 提取后验结果
mpost = convertToCodaObject(models[[1]])

# 查看系统发育信号
summary(mpost$Rho)$quant

mpost = convertToCodaObject(models[[2]])

summary(mpost$Rho)$quant