##################################################################################################
# 1. Set the base directory using your favorite method
# 1. 使用你喜欢的方式设置基础工作目录
##################################################################################################
# setwd("...")


##################################################################################################
# 2. INPUT AND OUTPUT OF THIS SCRIPT (BEGINNING)
# 2. 本脚本的输入和输出（开始）
##################################################################################################
#   INPUT. Original datafiles of the case study, placed in the data folder.
#   2.1 输入：研究案例的原始数据文件，放在 data 文件夹中。

#   OUTPUT. Unfitted models, i.e., the list of Hmsc model(s) that have been defined
#   2.2 输出：未拟合的模型（已定义但尚未拟合的 Hmsc 模型列表）

#           but not fitted yet, stored in the file "models/unfitted_models.RData".
#           2.3 这些模型将存储在文件 models/unfitted_models.RData 中。
##################################################################################################
# INPUT AND OUTPUT OF THIS SCRIPT (END)
# 本脚本的输入和输出（结束）
##################################################################################################


##################################################################################################
# 3. MAKE THE SCRIPT REPRODUCIBLE (BEGINNING)
# 3. 让脚本具有可重复性（开始）
##################################################################################################
set.seed(1)
# 3.1 设置随机数种子为 1（保证每次运行结果一致）
##################################################################################################
# MAKE THE SCRIPT REPRODUCIBLE (END)
# 让脚本具有可重复性（结束）
##################################################################################################


##################################################################################################
# 4. LOAD PACKAGES (BEGINNING)
# 4. 加载 R 包（开始）
##################################################################################################
library(Hmsc)
# 4.1 加载 Hmsc 包（用于构建分层多物种模型）
##################################################################################################
# LOAD PACKAGES (END)
# 4. 加载 R 包（结束）
##################################################################################################


##################################################################################################
# 5. SET DIRECTORIES (BEGINNING)
# 5. 设置目录（开始）
##################################################################################################
localDir = "."
# 5.1 将当前目录设置为项目目录

dataDir = file.path(localDir, "data")
# 5.2 定义数据文件夹路径 data/

modelDir = file.path(localDir, "models")
# 5.3 定义模型文件夹路径 models/

if(!dir.exists(modelDir)) dir.create(modelDir)
# 5.4 如果 models 文件夹不存在，则创建它
##################################################################################################
# SET DIRECTORIES (END)
# 5. 设置目录（结束）
##################################################################################################


##################################################################################################
# 6. READ AND EXPLORE THE DATA (BEGINNING)
# 6. 读取并探索数据（开始）
##################################################################################################
# Write here the code needed to read in the data, and explore it by (with View, plot, hist, ...)
# 6.1 在此写入加载数据的代码，并使用 View、plot、hist 等函数查看数据结构

# to get and idea of it and ensure that the data are consistent
# 6.2 以便你理解数据内容并确保数据一致且无错误
##################################################################################################
# READ AND EXPLORE THE DATA (END)
# 6. 读取并探索数据（结束)
##################################################################################################


##################################################################################################
# 7. SET UP THE MODEL (BEGINNING)
# 7. 设置 HMSC 模型（开始）
##################################################################################################
# Note that many of the components are optional 
# 7.1 注意：很多模型组件是可选的

# Here instructions are given at generic level, the details will depend on the model
# 7.2 下列说明是通用指导，具体内容依你的模型需求而定


# Organize the community data in the matrix Y
# 7.3 将群落数据整理为矩阵 Y（样地 × 物种）


# Organize the environmental data into a dataframe XData
# 7.4 将环境变量整理为数据框 XData

# Define the environmental model through XFormula
# 7.5 使用 XFormula 定义环境变量公式
# XFormula = ~ ...


# Organize the trait data into a dataframe TrData
# 7.6 将性状数据整理为数据框 TrData

# Define the trait model through TrFormula
# 7.7 使用 TrFormula 定义性状模型
# TrFormula = ~ ..


# Set up a phylogenetic (or taxonomic tree) as myTree
# 7.8 设置系统发育树 or 分类树，命名为 myTree


# Define the studyDesign as a dataframe 
# 7.9 将 studyDesign 定义为数据框（包含随机效应分组）

# For example, if you have sampled the same locations over multiple years, you may define
# 7.10 举例：若在多个年份重复采样同一地点，可以这样定义：

# studyDesign = data.frame(sample = ..., year = ..., location = ...)
# studyDesign = data.frame(sample = ..., year = ..., location = ...)


# Set up the random effects
# 7.11 设置随机效应（random effects）



# For example, you may define year as an unstructured random effect
# 7.12 例如，可将 year 定义为非结构化随机效应

# rL.year = HmscRandomLevel(units = levels(studyDesign$year))
# rL.year = HmscRandomLevel(units = levels(studyDesign$year))



# For another example, you may define location as a spatial random effect
# 7.13 另一个例子：将 location 定义为空间随机效应

# rL.location = HmscRandomLevel(sData = locations.xy)
# rL.location = HmscRandomLevel(sData = locations.xy)

# Here locations.xy would be a matrix (one row per unique location)
# 7.14 其中 locations.xy 是矩阵（每行表示一个地点）

# where row names are the levels of studyDesign$location,
# 7.15 行名必须对应 studyDesign$location 的水平

# and the columns are the xy-coordinates
# 7.16 列为该地点的 XY 坐标



# For another example, you may define the sample = sampling unit = row of matrix Y
# 7.17 另一例：可将 sample（等于 Y 的行）作为随机效应

# as a random effect, in case you are interested in co-occurrences at that level
# 7.18 若你想研究样地层级的物种共现，则可这样设定

# rL.sample = HmscRandomLevel(units = levels(studyDesign$sample))
# rL.sample = HmscRandomLevel(units = levels(studyDesign$sample))



# Use the Hmsc model constructor to define a model
# 7.19 使用 Hmsc 构造函数创建模型

# m = Hmsc(Y=Y,
# 7.20 构建模型的主函数 Hmsc

#          distr="probit",
# 7.21 设定分布为 probit（适用于 presence–absence 数据）

#          XData = XData,  XFormula=XFormula,
# 7.22 输入环境变量及其公式

#          TrData = TrData, TrFormula = TrFormula,
# 7.23 输入性状变量及其公式

#          phyloTree = myTree,
# 7.24 输入系统发育树

#          studyDesign = studyDesign, 
# 7.25 输入研究设计（样地、年份、地点等分组）

#          ranLevels=list(year=rL.year, location = rL.location, sample = rL.sample))
# 7.26 输入所有随机效应（year, location, sample）



# note that in the random effects the left-hand sides in the list (year, location, sample)
# 7.27 注意：随机效应列表中的名称（year, location, sample）

# refer to the columns of the studyDesign
# 7.28 必须与 studyDesign 中的列名一致



# In this example we assumed the probit distribution as appropriate for presence-absence data
# 7.29 此处使用 probit 分布，因为适用于 presence–absence 数据


# It is always a good idea to look at the model object, so type m to the console and press enter
# 7.30 建议你查看模型对象：在控制台输入 m 并按回车

# Look at the components of the model by exploring m$...
# 7.31 使用 m$... 查看模型内部组件
##################################################################################################
# SET UP THE MODEL (END)
# 7. 设置 HMSC 模型（结束）
##################################################################################################


##################################################################################################
# 8. COMBINING AND SAVING MODELS (START)
# 8. 合并并保存模型（开始）
##################################################################################################
# models = list(m, m.alternative)
# 8.1 将主模型和备选模型放入列表 models = list(m, m.alternative)

# names(models) = c("my.main.model","my.alternative model")
# 8.2 给每个模型命名

# save(models, file = file.path(modelDir, "unfitted_models.RData"))
# 8.3 将未拟合的模型保存到 models/unfitted_models.RData
##################################################################################################
# COMBINING AND SAVING MODELS (END)
# 8. 合并并保存模型（结束）
##################################################################################################


##################################################################################################
# 9. TESTING THAT MODELS FIT WITHOUT ERRORS (START)
# 9. 测试模型是否能正常拟合（开始）
##################################################################################################
#for(i in 1:length(models)){
#  9.1 循环遍历每个模型
#  print(i)
#  9.2 打印当前模型编号
#  sampleMcmc(models[[i]],samples=2)
#  9.3 运行一个极小的 MCMC（仅 2 个样本）测试是否报错
#}
##################################################################################################
# TESTING THAT MODELS FIT WITHOUT ERRORS (END)
# 9. 测试模型是否能正常拟合（结束）
##################################################################################################
