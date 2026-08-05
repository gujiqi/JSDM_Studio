# JSDM Studio 参数详细解释手册

本手册对应 JSDM Studio 当前界面，逐项解释每一个常用参数、每个输入文件的行列含义、什么时候该选、什么时候要小心。

当前版本的建模核心是 **Hmsc R 包**。JSDM Studio 不修改 Hmsc 的统计模型，只是把常见 Hmsc/JSDM 工作流做成图形化、可重复的流程。

---

# 0. JSDM Studio 在做什么？

JSDM Studio 用于联合物种分布模型（Joint Species Distribution Models, JSDMs）。它处理的是多物种群落数据，而不是单物种数据。

最核心的数据结构是：

```text
Y      = 样方 × 物种 的响应矩阵
XData  = 样方 × 环境变量 的解释变量矩阵
TrData = 物种 × 性状 的性状矩阵，可选
studyDesign = 样方 × 分组变量 的随机效应设计表，可选
coordinates = 样方 × 坐标 的空间位置表，可选
phyloTree = 物种系统发育树，可选
```

模型主要回答：

```text
1. 哪些环境变量影响群落组成？
2. 不同物种对环境变量的响应是否不同？
3. 物种性状能否解释这些响应差异？
4. 控制环境后，物种之间是否还有残差关联？
5. 空间、地点、样地、年份等随机效应有多重要？
6. 模型拟合和预测能力如何？
7. 结果是否可重复、能否写进论文？
```

---

# 1. 数据上传页面

## 1.1 Y.csv：物种响应矩阵

### 含义

`Y.csv` 是最重要的文件，表示每个样方中每个物种的响应值。

示例：

| sample | sp1 | sp2 | sp3 |
|---|---:|---:|---:|
| site1 | 1 | 0 | 1 |
| site2 | 0 | 1 | 1 |
| site3 | 1 | 0 | 0 |

### 行是什么意思？

```text
每一行 = 一个样方 / 样点 / 样本 / 调查单元
```

例如：

```text
一个 1 m² 样方
一个调查点
一个样带上的一个样点
一个陷阱样本
一个时间-地点组合
```

### 列是什么意思？

```text
每一列 = 一个物种
```

例如：

```text
sp1
Sphagnum_palustre
Polytrichum_commune
species_001
```

### 单元格是什么意思？

取决于你选择的 `distr`：

| Y 里面的值 | 推荐 distr | 含义 |
|---|---|---|
| 0/1 | probit | 物种不出现/出现 |
| 0,1,2,3... | poisson | 个体数、记录次数、株数等计数 |
| 连续数值 | normal | 盖度、生物量、连续丰度指数等 |

### 注意

1. 第一列一般作为行名，也就是样方 ID。
2. 行数必须和 `XData.csv` 一致。
3. 物种列名最好不要有空格、括号、中文标点。
4. 如果用 `TrData`，`Y` 的列名必须和 `TrData` 的行名对应。
5. 如果用系统发育树，`Y` 的列名必须和树的 tip labels 对应。
6. 不能有完全空的物种列。
7. 过多稀有物种会影响模型稳定性，正式分析前可考虑稀有种筛选。

---

## 1.2 XData.csv：环境变量矩阵

### 含义

`XData.csv` 是每个样方的环境变量、处理变量或解释变量。

示例：

| sample | pH | moisture | canopy | elevation | substrate |
|---|---:|---:|---:|---:|---|
| site1 | 5.2 | 0.73 | 65 | 1200 | soil |
| site2 | 4.8 | 0.81 | 72 | 1280 | log |
| site3 | 6.1 | 0.42 | 30 | 900 | rock |

### 行是什么意思？

```text
每一行 = 一个样方
```

必须与 `Y.csv` 的样方顺序或样方 ID 对应。

### 列是什么意思？

```text
每一列 = 一个环境变量
```

常见变量：

```text
pH
moisture
canopy
elevation
temperature
precipitation
substrate
habitat
soil_type
land_use
```

### 数值变量

例如：

```text
pH
moisture
canopy
elevation
temperature
```

可以直接进入公式：

```r
~ pH + moisture + canopy + elevation
```

### 分类变量

例如：

```text
substrate = soil / log / rock
habitat = forest / grassland / wetland
```

R/Hmsc 会把它们作为 factor 处理。

### 注意

1. 缺失值可能导致模型失败。
2. 分类变量的拼写要统一，例如 `soil` 和 `Soil` 会被当成不同水平。
3. 连续变量量纲差异很大时，建议标准化。
4. 变量之间高度相关时，参数解释要谨慎。
5. 变量名不要带空格，建议用下划线。

---

## 1.3 TrData.csv：物种性状矩阵，可选

### 含义

`TrData.csv` 是物种层面的性状数据。

示例：

| species | height_mm | life_form | reproduction |
|---|---:|---|---|
| sp1 | 20 | mat | sexual |
| sp2 | 55 | turf | asexual |
| sp3 | 12 | cushion | sexual |

### 行是什么意思？

```text
每一行 = 一个物种
```

行名必须对应 `Y.csv` 的物种列名。

### 列是什么意思？

```text
每一列 = 一个性状变量
```

例如：

```text
height
seed_mass
life_form
body_size
feeding_guild
dispersal_type
reproductive_mode
```

### 什么时候用？

当你想回答：

```text
物种性状是否能解释物种对环境的响应差异？
哪些性状对应更强的环境响应？
```

### 注意

1. 物种名必须和 `Y` 的列名对应。
2. 性状缺失值需要提前处理。
3. 性状太多、物种太少时，模型可能不稳定。
4. 分类性状会产生多个水平参数。
5. 初学者第一次测试可以不上传 TrData。

---

## 1.4 studyDesign.csv：随机效应设计表，可选

### 含义

`studyDesign.csv` 表示每个样方属于哪个分组，用于随机效应。

示例：

| sample | site | plot | year |
|---|---|---|---|
| s1 | A | A1 | 2024 |
| s2 | A | A2 | 2024 |
| s3 | B | B1 | 2025 |

### 行是什么意思？

```text
每一行 = 一个样方
```

### 列是什么意思？

```text
每一列 = 一个分组变量
```

例如：

```text
site
plot
transect
region
year
observer
sample
```

### 什么时候用？

如果数据不是完全独立，应该考虑随机效应：

```text
多个样方在同一个地点
多个样方在同一条样带
同一地点多年重复调查
同一实验区有多个样方
```

### 注意

1. 如果你在界面里选择随机效应列名为 `site`，那么 `studyDesign.csv` 中必须有 `site` 列。
2. 如果没有上传 studyDesign，程序可能会自动构造简单的 sample-level 随机效应。
3. 随机效应太多、样本太少会导致模型慢或不稳定。

---

## 1.5 coordinates.csv：空间坐标，可选

### 含义

`coordinates.csv` 表示样方的空间位置。

示例：

| sample | longitude | latitude |
|---|---:|---:|
| s1 | 120.123 | 31.123 |
| s2 | 120.135 | 31.130 |

或：

| sample | x | y |
|---|---:|---:|
| s1 | 10.2 | 5.4 |
| s2 | 12.1 | 6.2 |

### 行是什么意思？

```text
每一行 = 一个样方
```

### 列是什么意思？

```text
经度/纬度 或 平面坐标 x/y
```

### 什么时候用？

如果你认为样方之间存在空间自相关：

```text
距离近的样方更相似
物种分布有空间结构
环境变量不能完全解释空间聚集
```

### 注意

1. 坐标必须和 `Y` 的样方对应。
2. 经纬度和投影坐标不要混用。
3. 空间模型会更慢。
4. 小数据第一次测试不建议先开空间随机效应。

---

## 1.6 phyloTree：系统发育树，可选

### 含义

系统发育树表示物种之间的进化关系。

### 什么时候用？

如果你想回答：

```text
亲缘关系近的物种是否有相似环境响应？
环境响应是否存在系统发育信号？
```

### 注意

1. 树的 tip labels 必须和 `Y` 的物种列名对应。
2. 树格式通常是 Newick。
3. 物种名不匹配是最常见错误。
4. 不确定时第一次不要上传树。

---

# 2. 模型参数页面

## 2.1 distr：响应变量分布

### 含义

`distr` 告诉模型 Y 中的数据是什么类型。

### 可选项

| distr | 适合数据 | 例子 |
|---|---|---|
| probit | 0/1 出现-不出现 | presence/absence |
| poisson | 非负整数计数 | 个体数、株数、记录次数 |
| normal | 连续数据 | 盖度、生物量、连续丰度 |

### 怎么选？

```text
Y 只有 0 和 1 → probit
Y 是 0,1,2,3... → poisson
Y 是连续数值 → normal
```

### 注意

1. `probit` 不适合计数。
2. `poisson` 不适合负数或小数。
3. `normal` 不适合纯 0/1 数据。
4. 盖度百分比有时不严格符合 normal，解释要谨慎。

---

## 2.2 XFormula：环境变量公式

### 含义

`XFormula` 决定哪些环境变量进入模型。

### 常见写法

使用所有 XData 变量：

```r
~ .
```

只使用部分变量：

```r
~ pH + moisture + canopy + elevation
```

加入二次项：

```r
~ pH + I(pH^2) + moisture
```

加入交互项：

```r
~ pH * moisture
```

### 行列含义

公式里的变量名必须来自 `XData.csv` 的列名。

### 注意

1. 变量名必须拼写完全一致。
2. `~ .` 很方便，但会把 XData 中所有列都放进去，可能包括不该进入模型的 ID 列。
3. 变量过多、样方过少会导致模型不稳定。
4. 强相关变量同时进入模型时，参数解释会变困难。
5. 测试时建议先用少数关键变量。

---

## 2.3 使用物种性状 TrData

### 含义

是否把物种性状放进模型。

勾选后，模型可以估计：

```text
性状如何影响物种对环境变量的响应
```

### 什么时候勾选？

如果你上传了 `TrData.csv`，并且想解释物种响应差异。

### 不建议勾选的情况

```text
TrData 缺失很多
物种数很少
性状变量太多
只是想先测试软件是否能跑通
```

---

## 2.4 TrFormula：性状公式

### 含义

决定哪些性状进入模型。

使用所有性状：

```r
~ .
```

使用部分性状：

```r
~ height_mm + life_form + reproduction
```

### 注意

1. 变量名来自 `TrData.csv` 的列名。
2. 分类性状会展开成多个参数。
3. 性状和环境响应之间的关系通常通过 Gamma 参数解释。

---

## 2.5 使用 phylogeny

### 含义

是否加入系统发育结构。

### 什么时候用？

如果你上传了系统发育树，且有明确问题：

```text
物种响应是否受亲缘关系影响？
亲缘关系近的物种是否更相似？
```

### 注意

1. 树 tip labels 必须与物种名匹配。
2. 系统发育模型会增加复杂度。
3. 不确定时先不要勾选。

---

# 3. 随机效应与空间参数

## 3.1 random_mode：随机效应模式

### 常见选项

| 模式 | 含义 |
|---|---|
| none | 不使用随机效应 |
| sample | 使用样方/分组随机效应 |
| spatial | 使用空间随机效应 |

### 怎么选？

```text
第一次测试 → sample 或 none
有空间坐标且关注空间自相关 → spatial
数据很简单 → none
```

### 注意

随机效应用于处理非独立性，但会增加模型复杂度和运行时间。

---

## 3.2 random_effect_column：随机效应列名

### 含义

告诉程序使用 studyDesign 中哪一列作为随机效应。

例如：

```text
site
plot
region
sample
year
```

### 例子

如果 `studyDesign.csv` 是：

| sample | site | year |
|---|---|---|
| s1 | A | 2024 |
| s2 | A | 2024 |
| s3 | B | 2025 |

你填：

```text
site
```

表示同一个 site 内的样方共享一个随机效应。

### 注意

1. 列名必须存在。
2. 如果没有上传 studyDesign，填 `sample` 通常表示样方级随机效应。
3. 分组水平太多或太少都可能影响模型。

---

## 3.3 longitude_column / latitude_column：坐标列名

### 含义

告诉程序 coordinates 文件中哪两列是空间坐标。

常见：

```text
longitude / latitude
x / y
lon / lat
```

### 注意

1. 列名必须和 coordinates.csv 完全一致。
2. 经度纬度适合大范围空间位置，但距离计算要谨慎。
3. 平面投影坐标更适合空间距离建模。

---

## 3.4 spatial_method：空间方法

### 常见选项

| 方法 | 含义 |
|---|---|
| Full | 完整空间随机效应 |
| NNGP | 近邻高斯过程，适合较大数据 |
| GPP | 预测过程近似，适合较大数据 |

### 怎么选？

```text
小数据测试 → Full
样方很多 → NNGP
不熟悉空间模型 → 先不要 spatial
```

### 注意

空间方法会显著增加运行时间。第一次跑通软件时不建议先开空间模型。

---

## 3.5 nNeighbours：NNGP 邻居数

### 含义

当使用 NNGP 空间近似时，每个样方参考多少个最近邻居。

### 常见值

```text
5
10
15
20
```

### 怎么选？

```text
测试 → 10
样方较多 → 10 到 20
```

### 注意

1. 值越大，近似越精细，但越慢。
2. 只在 NNGP 等近似空间方法中重要。
3. 使用 Full 时通常不用关心。

---

# 4. MCMC 参数

## 4.1 samples：保留样本数

### 含义

MCMC 结束后保留多少个后验样本用于推断。

### 测试值

```text
20
50
100
```

### 正式分析

通常需要更大，例如：

```text
1000
2000
5000
```

具体要看模型复杂度和收敛诊断。

### 注意

1. `samples` 越大，结果越稳定，但运行越慢。
2. 测试小 samples 只能证明程序能跑，不能作为论文结果。
3. 正式分析必须查看收敛。

---

## 4.2 transient：burn-in / 预热期

### 含义

MCMC 前期丢弃的迭代数，不用于最终结果。

### 测试值

```text
10
50
100
```

### 正式分析

可能需要：

```text
1000
5000
```

或更多，取决于收敛情况。

### 注意

1. transient 太小可能保留未收敛阶段。
2. transient 不是越大越好，太大增加时间。
3. 应结合 traceplot 和 ESS 判断。

---

## 4.3 thin：抽稀间隔

### 含义

每隔多少次迭代保存一个样本。

例如：

```text
thin = 1：每次都保存
thin = 10：每 10 次保存一次
```

### 注意

1. thin 越大，总迭代通常越多。
2. thin 可减少样本自相关和文件大小。
3. 测试时用 1 即可。

---

## 4.4 nChains：链数

### 含义

MCMC 独立链的数量。

### 测试值

```text
2
```

### 正式分析建议

```text
至少 2
常用 4
```

### 注意

1. 多链可以检查 Gelman/Rhat/PSRF 收敛。
2. 只有 1 条链时，很多收敛诊断不可靠。
3. 链数越多，计算量越大。

---

## 4.5 nParallel：并行数

### 含义

同时运行多少个链或任务。

### 建议

```text
新手 Windows → 1
电脑较强 → 2 或 4
```

### 注意

1. 不要超过 CPU 核心数。
2. Windows 下并行有时更容易出问题。
3. 如果报错，先改回 1。

---

## 4.6 verbose：进度显示间隔

### 含义

Hmsc 在控制台中显示 MCMC 进度的间隔。

### 例子

```text
verbose = 5：更频繁显示
verbose = 50：较少显示
verbose = 0：可能不显示
```

### 注意

GUI 进度条无法实时显示 sampleMcmc 内部每一步，控制台 verbose 更接近底层进度。

---

## 4.7 seed：随机种子

### 含义

控制随机数，使结果更可重复。

### 建议

```text
123
2024
1
```

### 注意

1. 同一数据、同一参数、同一 seed，结果更容易复现。
2. 不同 seed 可用于检查结果稳定性。
3. 论文分析应记录 seed。

---

# 5. 输出设置参数

## 5.1 保存模型

### 输出

```text
models/hmsc_model.rds
models/models_thin_*_samples_*_chains_*.RData
```

### 含义

保存拟合后的模型对象，便于以后重新读取、继续分析、复现结果。

### 建议

正式分析一定勾选。

---

## 5.2 computePredictedValues：计算预测值

### 输出

```text
results/predicted_values.rds
```

### 含义

计算模型对训练数据的预测值。

### 用途

后续模型拟合、预测能力、图表都可能需要它。

---

## 5.3 evaluateModelFit：解释性模型拟合

### 输出

```text
results/model_fit_explanatory.txt
tables/model_fit_explanatory_*.csv
plots/model_fit_explanatory_vs_predictive.pdf
```

### 含义

评估模型对已有数据的解释能力。

### 注意

解释性拟合不是独立预测能力，通常会比交叉验证更乐观。

---

## 5.4 交叉验证预测能力

### 输出

```text
results/model_fit_cross_validation.txt
tables/model_fit_cross_validation_*.csv
```

### 含义

通过把数据分成训练集/测试集，评估模型预测未见数据的能力。

### 注意

1. 比解释性拟合更接近真实预测能力。
2. 更耗时。
3. 小数据下结果波动较大。

---

## 5.5 nfolds：交叉验证折数

### 含义

交叉验证分成几折。

### 常见值

```text
2：最快，适合测试
5：常用
10：更常见但更慢
```

### 注意

nfolds 越大，运行时间越长。

---

## 5.6 WAIC

### 含义

WAIC 是模型比较指标，用于比较不同模型的相对支持度。

### 注意

1. 单个模型的 WAIC 意义有限。
2. 更适合多模型比较。
3. 计算可能较慢。
4. 测试阶段可先关闭。

---

## 5.7 MCMC 收敛诊断

### 输出

```text
results/MCMC_convergence.txt
plots/MCMC_traceplots.pdf
```

### 含义

检查 MCMC 是否稳定、是否收敛。

### 重点看

```text
traceplot 是否像稳定毛毛虫
PSRF/Rhat 是否接近 1
ESS 是否足够
```

### 注意

正式论文必须检查收敛，不能只看模型跑完。

---

## 5.8 参数估计 Beta/Gamma/Omega

### Beta

环境变量对物种的影响。

```text
行/列通常涉及环境变量和物种
数值表示后验均值、支持度等
```

### Gamma

性状如何解释物种环境响应。

```text
行/列通常涉及性状和环境响应关系
```

### Omega

物种残差关联。

```text
行 = 物种
列 = 物种
数值 = 控制环境和随机效应后的残差关联
```

### 注意

Omega 不等于直接种间竞争。它只是残差相关，可能来自未测环境变量、空间结构、相互作用或采样过程。

---

## 5.9 方差分解

### 输出

```text
tables/variance_partitioning_*.csv
plots/variance_partitioning.pdf
```

### 含义

展示不同解释组对物种变异的贡献。

### 解释

```text
某个物种的变异有多少由环境变量解释？
有多少由随机效应或空间结构解释？
```

---

## 5.10 物种关联 Omega

### 输出

```text
results/species_associations.rds
tables/Omega_*_*.csv
plots/Omega_associations_*.pdf
```

### 行列含义

```text
行 = 物种 A
列 = 物种 B
单元格 = 两个物种的残差关联强度
```

### 解释

正值：

```text
两个物种在控制环境后仍倾向共同出现
```

负值：

```text
两个物种在控制环境后倾向不共同出现
```

### 注意

不能直接写成竞争或互利，除非有额外实验或生态证据。

---

## 5.11 环境梯度预测图

### 输出

```text
plots/predictions_environmental_gradients.pdf
```

### 含义

展示沿某个环境变量变化时，物种预测出现概率或响应值如何变化。

### 注意

1. 连续变量更适合画梯度。
2. 分类变量解释要谨慎。
3. 预测图依赖模型设定和变量范围。
4. 不要外推到数据范围之外。

---

## 5.12 HTML 报告

### 输出

```text
HmscGUI_report.html
```

或重命名后可能显示：

```text
JSDM_Studio_report.html
```

### 含义

自动生成结果索引和简要说明，便于查看输出文件。

---

# 6. 运行页面参数

## 检查数据

### 含义

读取上传的文件，检查行列维度、基本格式和是否能进入后续流程。

### 注意

通过检查不代表统计模型一定合理，只代表基本格式可读。

---

## 运行 Hmsc 并生成结果

### 含义

根据当前界面参数：

```text
1. 保存 used_config.yml
2. 构建 Hmsc 模型
3. 运行 MCMC
4. 计算勾选的输出
5. 保存表格、图、模型对象和报告
```

---

## 下载结果 ZIP

### 含义

把当前运行输出目录打包，便于保存、分享或论文复现。

---

# 7. 测试参数 vs 正式论文参数

## 测试参数

```text
samples = 20
transient = 10
thin = 1
nChains = 2
nParallel = 1
```

用途：

```text
检查软件能否跑通
检查数据格式
检查输出文件是否生成
```

不能用于正式论文。

## 正式分析参数

没有固定唯一值，要根据收敛诊断决定。通常需要：

```text
更大的 samples
更大的 transient
至少 2 到 4 条链
检查 traceplot
检查 ESS / PSRF / Rhat
```

论文中应报告：

```text
distr
XFormula
TrFormula
随机效应结构
MCMC samples/transient/thin/nChains
seed
模型诊断结果
预测或交叉验证结果
```
