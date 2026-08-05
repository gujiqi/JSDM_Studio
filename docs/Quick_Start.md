# JSDM Studio 快速开始

## 普通用户

1. 安装 R for Windows  
   https://cran.r-project.org/bin/windows/base/

2. 安装 `JSDMStudio_Setup.exe`

3. 双击桌面图标：

```text
JSDM Studio
```

4. 第一次启动请等待 R 包自动检查/安装。

5. 打开界面后，先用示例数据测试。

---

## 第一次测试参数

数据上传：

```text
Y.csv = examples/Y.csv
XData.csv = examples/XData.csv
traits = 不选
studyDesign = 不选
coordinates = 不选
phyloTree = 不选
```

模型参数：

```text
distr = probit
XFormula = ~ pH + moisture + canopy + elevation
random_mode = sample
random_effect_column = sample
```

MCMC 参数：

```text
samples = 20
transient = 10
thin = 1
nChains = 2
nParallel = 1
verbose = 5
seed = 123
```

第一次输出只勾选：

```text
保存模型
计算预测值
解释性模型拟合
MCMC 收敛诊断
生成 HTML 报告
```

先不要勾选：

```text
WAIC
环境梯度预测
物种关联 Omega
```

基础跑通后再逐步打开高级输出。

---

## 结果在哪里？

运行成功后，结果在：

```text
output/hmsc_日期时间/
```

主要文件：

```text
used_config.yml
RUN_COMPLETE.txt
JSDM Studio_report.html
models/
results/
tables/
plots/
```

---

## 出错怎么办？

优先看：

```text
startup_log.txt
```

如果独立窗口打不开，用备用模式：

```text
Launch_in_default_browser.bat
```
