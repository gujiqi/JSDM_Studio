# JSDM Studio / JSDM Studio：图标美化 + 结果增强版

这是一个基于 Hmsc R 包的本地图形化工作台。用户需要安装 R，但不需要写 R 代码。

## 使用方式

第一次：

1. 双击 `install_packages.bat`
2. 等待 R 包安装完成

以后：

1. 双击 `run_app.bat`
2. 浏览器会自动打开图形界面
3. 上传数据
4. 设置模型参数
5. 点击“运行 Hmsc 并生成结果”
6. 下载结果 ZIP

## 第一次测试推荐参数

- Y：`Y.csv`
- XData：`XData.csv`
- distr：`probit`
- XFormula：`~ pH + moisture + canopy`
- samples：20
- transient：10
- thin：1
- nChains：2
- nParallel：1

测试阶段可以关闭 WAIC。交叉验证、环境梯度预测和物种关联可能比较耗时。

## 新增结果输出

本版本参考 Hmsc 官方脚本式工作流，把结果部分扩展为多个文件夹：

```text
output/hmsc_日期时间/
├── used_config.yml
├── RUN_COMPLETE.txt
├── JSDM Studio_report.html
├── inputs/
│   └── 本次上传的原始输入文件副本
├── models/
│   ├── unfitted_models.RData
│   ├── unfitted_model.rds
│   ├── models_thin_*_samples_*_chains_*.RData
│   └── hmsc_model.rds
├── results/
│   ├── model_structure.txt
│   ├── predicted_values.rds
│   ├── model_fit_explanatory.txt
│   ├── model_fit_cross_validation.txt
│   ├── WAIC.txt
│   ├── MCMC_convergence.txt
│   ├── parameter_estimates.txt
│   └── species_associations.rds
├── tables/
│   ├── data_dimensions.csv
│   ├── species_summary.csv
│   ├── environment_summary.csv
│   ├── model_fit_explanatory_*.csv
│   ├── model_fit_cross_validation_*.csv
│   ├── parameter_estimates_Beta_*.csv
│   ├── parameter_estimates_Gamma_*.csv
│   ├── variance_partitioning_*.csv
│   └── Omega_*_*.csv
└── plots/
    ├── MCMC_traceplots.pdf
    ├── model_fit_explanatory_vs_predictive.pdf
    ├── Beta_plot.pdf
    ├── Gamma_plot.pdf
    ├── variance_partitioning.pdf
    ├── Omega_associations_*.pdf
    └── predictions_environmental_gradients.pdf
```

## 输出解释

- `used_config.yml`：本次运行使用的全部参数，保证可重复。
- `models/unfitted_models.RData`：未拟合模型，类似 S1_define_models_template.R 的输出。
- `models/models_thin_*_samples_*_chains_*.RData`：拟合后的模型，类似 S2_fit_models.R 的命名风格。
- `results/MCMC_convergence.txt` 和 `plots/MCMC_traceplots.pdf`：收敛诊断。
- `results/model_fit_explanatory.txt`：解释性模型拟合。
- `results/model_fit_cross_validation.txt`：交叉验证预测能力。
- `tables/parameter_estimates_Beta_*.csv`：环境响应参数 Beta。
- `tables/parameter_estimates_Gamma_*.csv`：性状影响参数 Gamma。
- `plots/variance_partitioning.pdf`：方差分解图。
- `plots/Omega_associations_*.pdf`：物种残差关联图。
- `plots/predictions_environmental_gradients.pdf`：环境梯度预测图。
- `JSDM Studio_report.html`：自动生成的结果索引报告。

## 注意

正式论文分析不要使用测试小参数。需要提高 `samples`、`transient`，并检查 MCMC 诊断。

## 进度与剩余时间显示

本版本新增“进度与预计剩余时间”面板：

- 运行前会根据样方数、物种数、MCMC 参数和输出选项估算总耗时。
- 运行时显示当前阶段、进度百分比、已用时间和预计剩余时间。
- MCMC 的 `sampleMcmc()` 是一个长时间阻塞步骤，GUI 不能像控制台 verbose 那样逐次刷新内部迭代，只能在进入和完成 MCMC 阶段时更新。控制台仍会显示 Hmsc 自己的 verbose 信息。
- 预计剩余时间是粗略估计，正式大数据会受电脑性能、随机效应、空间模型、交叉验证和输出图影响。

## 数字参数滑块 + 自定义输入框

本版本把常用数字参数改成“左侧滑块 + 右侧数字输入框”：

- 新手可以拖动滑块；
- 高级用户可以在右侧直接输入精确自定义数值；
- 如果右侧填写超过滑块最大范围，运行时会使用右侧数字框的真实数值；
- 适用于 samples、transient、thin、nChains、nParallel、verbose、seed、nfolds、support 阈值等参数。


## 用户安装说明

详细说明见：

```text
docs/用户安装说明_必读.md
docs/快速开始.md
用户安装最短说明.txt
```

当前版本不是“只装 exe 就完全够”的版本。用户需要先安装 R for Windows；WebView2 Runtime 大多数 Windows 10/11 已自带；R 包由程序自动检查和安装。


## 如何从 ZIP 生成 Setup.exe

详细步骤见：

```text
从ZIP生成Setup安装包_详细说明.md
生成安装包检查表.txt
docs/从ZIP生成Setup安装包_详细说明.md
docs/生成安装包检查表.txt
```

最短顺序：

```text
1. 进入 JSDMStudio
2. 双击 install_packages.bat
3. 双击 Build_WebView2_Launcher.bat
4. 双击 Start_WebView2_Window.bat 测试
5. 进入 JSDMStudio\installer
6. 双击 build_installer.bat
7. 到 installer\output 找 JSDMStudio_Setup.exe
```


## 参数与结果详细解释

本版本新增了三个详细手册：

```text
参数详细解释手册.md
运行结果文件详细解释.md
完整用户解释手册_参数与结果.md
```

同样文件也保存在：

```text
docs/参数详细解释手册.md
docs/运行结果文件详细解释.md
docs/完整用户解释手册_参数与结果.md
```

内容包括：

```text
每个输入文件的行列含义
每个模型参数怎么填
每个 MCMC 参数的意义
每个输出选项产生什么文件
results/tables/plots 中每个主要文件的解释
每个表格行列是什么意思
哪些结果可以写进论文
哪些解释需要谨慎
```
