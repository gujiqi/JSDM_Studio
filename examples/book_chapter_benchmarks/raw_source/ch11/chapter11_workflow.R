# =========================================================
# 第11章：芬兰鸟类案例（HMSC）
# =========================================================

# -------------------------
# 1. 设置工作目录
# -------------------------
wd = "C:/Users/Google/Documents/R/Joint Species Distribution Modelling in R/section_11_1_birds_2020_05_31/Section_11_1_birds"
setwd(wd)

localDir = "."
data.directory = file.path(localDir, "data")
model.directory = file.path(localDir, "models")

# -------------------------
# 2. 载入包并读入数据
# -------------------------
library(Hmsc)
library(ape)
set.seed(1)

da = read.csv(file.path(data.directory, "data.csv"), stringsAsFactors = TRUE)
da$Route = as.factor(da$Route)

# 只使用 2014 年数据
da = droplevels(subset(da, Year == 2014))

# 环境变量
XData = data.frame(
  Route = da$Route,
  hab = da$Habitat,
  clim = da$AprMay
)

# 群落矩阵：转为 presence-absence
Y = as.matrix(da[, -c(1:9)]) > 0
Y = apply(Y, MARGIN = 2, FUN = as.numeric)

# 空间坐标
xy = as.matrix(cbind(da$x, da$y))
rownames(xy) = da$Route
colnames(xy) = c("x-coordinate", "y-coordinate")

# 物种性状
alltraits = read.csv(file.path(data.directory, "traits.csv"), stringsAsFactors = TRUE)
Species = alltraits$Species
TrData = data.frame(
  Species = alltraits$Species,
  row.names = Species,
  Migration = alltraits$Migration,
  LogMass = log(alltraits$Mass)
)

# 系统发育树
phyloTree = read.tree(file.path(data.directory, "CTree.tre"))

# -------------------------
# 3. 看一下原始数据
# -------------------------
dim(Y)
head(Y[, c(1, 25, 50)])

S = rowSums(Y)
P = colMeans(Y)

head(XData)
head(xy)
head(TrData)
phyloTree

# =========================================================
# 11.1 Step 1：定义并拟合模型
# =========================================================

# studyDesign 中 route 要与 ranLevels 的名字一致
studyDesign = data.frame(route = XData$Route)

# 空间随机效应
rL = HmscRandomLevel(sData = xy)

# 固定效应公式
XFormula = ~ hab + poly(clim, degree = 2, raw = TRUE)

# 性状公式
TrFormula = ~ Migration + LogMass

# FULL 模型：环境 + 空间
m.FULL = Hmsc(
  Y = Y,
  XData = XData,
  XFormula = XFormula,
  phyloTree = phyloTree,
  TrData = TrData,
  TrFormula = TrFormula,
  distr = "probit",
  studyDesign = studyDesign,
  ranLevels = list(route = rL)
)

# ENV 模型：只有环境
m.ENV = Hmsc(
  Y = Y,
  XData = XData,
  XFormula = XFormula,
  phyloTree = phyloTree,
  TrData = TrData,
  TrFormula = TrFormula,
  distr = "probit"
)

# SPACE 模型：只有空间
m.SPACE = Hmsc(
  Y = Y,
  XData = XData,
  XFormula = ~1,
  phyloTree = phyloTree,
  TrData = TrData,
  TrFormula = TrFormula,
  distr = "probit",
  studyDesign = studyDesign,
  ranLevels = list(route = rL)
)

# MCMC 参数
nChains = 2
nParallel = 1
samples = 100
verbose = 100
transient = 50
thin = 1

models = list(m.FULL, m.ENV, m.SPACE)

# 拟合三个模型
for (i in 1:3) {
  models[[i]] = sampleMcmc(
    models[[i]],
    thin = thin,
    samples = samples,
    transient = transient,
    nChains = nChains,
    verbose = verbose,
    initPar = "fixed effects"
  )
}

# =========================================================
# 11.1 Step 2：检查 MCMC 收敛
# =========================================================

library(coda)

mpost = convertToCodaObject(
  models[[1]],
  spNamesNumbers = c(TRUE, FALSE),
  covNamesNumbers = c(TRUE, FALSE)
)

# Beta 参数的 Gelman-Rubin 诊断
psrf.beta = gelman.diag(mpost$Beta, multivariate = FALSE)$psrf

# Omega 参数抽样 200 个条目来检查
tmp = mpost$Omega[[1]]
z = dim(tmp[[1]])[2]
sel = sample(1:z, size = 200)

for (i in 1:length(tmp)) {
  tmp[[i]] = tmp[[i]][, sel]
}

psrf.omega = gelman.diag(tmp, multivariate = FALSE)$psrf

# =========================================================
# 11.1 Step 3：模型拟合优度与模型比较
# =========================================================

partition = createPartition(models[[1]], nfolds = 2, column = "route")

MF = list()
MFCV = list()

for (i in 1:3) {
  preds = computePredictedValues(models[[i]])
  MF[[i]] = evaluateModelFit(hM = models[[i]], predY = preds)
  
  preds = computePredictedValues(models[[i]], partition = partition)
  MFCV[[i]] = evaluateModelFit(hM = models[[i]], predY = preds)
}

# WAIC 比较
WAIC = unlist(lapply(models, FUN = computeWAIC))

# =========================================================
# 11.1 Step 4：参数解释
# =========================================================

# 看设计矩阵列名
head(models[[1]]$X)

# habitat 和 climate 两组做方差分解
groupnames = c("habitat", "climate")
group = c(1, 1, 1, 1, 1, 2, 2)

VP = list()
for (i in 1:2) {
  VP[[i]] = computeVariancePartitioning(
    models[[i]],
    group = group,
    groupnames = groupnames
  )
}

# Gamma 参数可视化
postGamma = getPostEstimate(models[[1]], parName = "Gamma")

par(mfrow = c(1, 2))
plotGamma(models[[1]], post = postGamma, param = "Support", supportLevel = 0.95)
plotGamma(models[[1]], post = postGamma, param = "Support", supportLevel = 0.85)

# 性状解释多少 Beta 变异
VP[[1]]$R2T$Beta

# 性状解释多少 Y 尺度变异
VP[[1]]$R2T$Y

# rho：系统发育信号
mpost = convertToCodaObject(models[[1]])
round(summary(mpost$Rho, quantiles = c(0.025, 0.5, 0.975))[[2]], 2)

# FULL 模型空间尺度参数 alpha
mpost = convertToCodaObject(models[[1]])
round(summary(mpost$Alpha[[1]], quantiles = c(0.025, 0.5, 0.975))[[2]][1:2, ], 2)

# SPACE 模型空间尺度参数 alpha
mpost = convertToCodaObject(models[[3]])
round(summary(mpost$Alpha[[1]], quantiles = c(0.025, 0.5, 0.975))[[2]][1:2, ], 2)

# =========================================================
# 11.1 Step 5：做预测
# =========================================================

grid = read.csv(file.path(data.directory, "grid_10000.csv"))
grid = droplevels(subset(grid, !(Habitat == "Ma")))

xy.grid = as.matrix(cbind(grid$x, grid$y))
XData.grid = data.frame(hab = grid$Habitat, clim = grid$AprMay)

m = models[[1]]
Gradient = prepareGradient(
  m,
  XDataNew = XData.grid,
  sDataNew = list(route = xy.grid)
)

predY = predict(m, Gradient = Gradient, predictEtaMean = TRUE)
EpredY = Reduce("+", predY) / length(predY)

# 物种 50：Corvus monedula
Cm = EpredY[, 50]

# 预测物种丰富度
S = rowSums(EpredY)

# 社区加权平均性状
CWM = (EpredY %*% m$Tr) / matrix(rep(S, m$nt), ncol = m$nt)

xy = grid[, 1:2]

# =========================================================
# 11.2 从后验中提取参数与不确定性
# =========================================================

m = models[[1]]

length(m$postList)
length(m$postList[[1]])

post = poolMcmcChains(m$postList)
length(post)

names(post[[1]])

# Beta 矩阵维度
dim(post[[1]]$Beta)

# X 矩阵列名
colnames(m$X)

# 第50个物种名称
m$spNames[50]

# 提取 beta[4, 50]
getvalue = function(p) {
  return(p$Beta[4, 50])
}

beta_4_50 = unlist(lapply(X = post, FUN = getvalue))

# 后验均值
mean(beta_4_50)

# 后验中位数
median(beta_4_50)

# 95% credible interval
quantile(beta_4_50, probs = c(0.025, 0.5, 0.975))

# Pr(beta > 0)
mean(beta_4_50 > 0)

# Pr(beta > 2)
mean(beta_4_50 > 2)

# 预测不确定性：以 species richness 为例
predY = predict(m, Gradient = Gradient, predictEtaMean = TRUE)
EpredY = Reduce("+", predY) / length(predY)

getS = function(p) {
  return(rowSums(p))
}

aS = simplify2array(lapply(X = predY, FUN = getS))
dim(aS)

ES = apply(aS, 1, mean)
sdS = sqrt(apply(aS, 1, var))
S25 = apply(aS > 25, 1, mean)

# =========================================================
# 11.4 用 HMSC 结果做 bioregionalisation
# =========================================================

RCP4 = kmeans(EpredY, 4)
RCP4$cluster = as.factor(RCP4$cluster)

RCP10 = kmeans(EpredY, 10)
RCP10$cluster = as.factor(RCP10$cluster)

# 选聚类数
library(factoextra)

p1 = fviz_nbclust(EpredY, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2) +
  labs(subtitle = "Elbow method")

p2 = fviz_nbclust(EpredY, kmeans, method = "silhouette") +
  labs(subtitle = "Silhouette method")

# 看 1 号和 43 号物种在各类中的中心值
RCP4$centers[, c(1, 43)]

# 每个区域最常见的物种
MCS = matrix(NA, m$ns, 4)

for (i in 1:4) {
  MCS[, i] = names(rev(sort(RCP4$centers[i, ])))
}

colnames(MCS) = c("Region 1", "Region 2", "Region 3", "Region 4")
MCS[1:5, ]

# =========================================================
# 11.5 与其他群落生态统计方法比较
# =========================================================

m = models[[1]]
Y = m$Y
XData = m$XData
TrData = m$TrData
xy = m$rL$route$s

# -------------------------
# 11.5.1 RDA 和方差分解
# -------------------------
library(vegan)

myrda = rda(Y ~ hab + clim, data = XData)
plot(myrda, display = c("bp", "site", "sp"))

library(ade4)
library(adespatial)

candidates = listw.candidates(xy, nb = c("gab"), weights = c("binary", "flin"))
modsel.Y = listw.select(
  Y,
  candidates,
  method = "FWD",
  MEM.autocor = "positive",
  p.adjust = TRUE
)

MEM.spe = modsel.Y$best$MEM.select

vY = vegdist(Y, method = "bray")
VP = varpart(vY, XData, MEM.spe)
plot(VP)

# -------------------------
# 11.5.2 Fourth-corner 分析
# -------------------------
library(ade4)

four2 = fourthcorner(XData, as.data.frame(Y), TrData, nrepet = 99, modeltype = 2)
summary(four2)

four4 = fourthcorner(XData, as.data.frame(Y), TrData, nrepet = 99, modeltype = 4)
summary(four4)

# -------------------------
# 11.5.3 共现分析
# -------------------------
library(EcoSimR)

summary(cooc_null_model(speciesData = t(Y)))

summary(cooc_null_model(speciesData = t(Y), algo = "sim2"))

# -------------------------
# 11.5.4 物种丰富度的单变量分析
# -------------------------
S = rowSums(Y)

m1 = glm(S ~ poly(clim, degree = 2) + hab, data = XData, family = "poisson")
m2 = glm(S ~ clim + hab, data = XData, family = "poisson")

AIC(m1, m2)

m3 = glm(S ~ poly(clim, degree = 2), data = XData, family = "poisson")
AIC(m1, m3)

plot(XData$clim, exp(predict(m3)))