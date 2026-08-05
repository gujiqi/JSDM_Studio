# JSDMWorkbench 全项目审计与修复报告

审计日期：2026-05-22

## 已修复

- 修复输出目录结构：每次运行都会创建 `inputs/`、`data/`、`models/`、`tables/`、`results/`、`plots/`、`predictions/`、`diagnostics/`、`workflow_scripts/`、`reproducible_script/`、`standard/`、`report/` 等目录。
- 修复下载逻辑：所有模型 ZIP 下载按钮现在都会返回真实 ZIP。未运行时生成 `not_run` 诊断 ZIP；检查失败时生成 `check_failed` 诊断 ZIP；不再写入 “No ZIP available” 文本假装 zip。
- 修复状态逻辑：scaffold-only 引擎不再显示 `Completed`，而是 `scaffold_only`；HMSC 真拟合失败显示 `fit_failed`；数据检查失败显示 `check_failed`。
- 修复 diagnostics：统一写入 `diagnostics/engine_status.json`、`diagnostics/data_check_messages.csv`、`diagnostics/session_info.txt`，错误时额外写入 `diagnostics/*_error.txt`。
- 修复 UI/server 对应：静态检查结果为 0 个缺失 UI input、0 个未使用 UI input、0 个缺失 output render/download。
- 修复 HMSC `hmsc_save_model`、jSDM `jsdm_out_config`、`jsdm_out_zip` 没有被后台读取的问题。
- 修复 HMSC `nfolds` 从错误层级读取的问题，现在读取 `outputs$predictions$nfolds`。
- 修复 `XData`/`TrData`/`studyDesign` 清洗：字符列会自动判断，数值型字符串转 numeric，分类字符串转 factor，逻辑列转 factor。
- 增强 Y 与 family 检查：binary/probit/binomial 必须 0/1；count/poisson 必须非负整数；beta/FC 必须在 0 和 1 之间；ordinal/OC/CAT 必须是正整数等级。
- 增强 formula 检查：HMSC、jSDM、GJAM、spOccupancy、sjSDM、boral 会检查公式是否合法，且公式变量必须存在于对应数据表。
- 修复 Compare Models：优先读取 `standard/run_summary.csv`、`standard/fit_metrics.csv`、`standard/effects_long.csv`、`standard/predictions_long.csv`、`standard/associations_long.csv`，不再只依赖可能不存在的非标准文件。
- 修复硬编码路径：移除 HMSC workflow templates 中作者电脑路径 `C:/Users/Google/Documents/R/...`。
- 修复 Windows 路径稳健性：`run_app.bat`、`install_packages.bat`、`run_cli.bat` 对 `%~dp0` 加引号，支持带空格路径。
- 修复 `install_packages.R` 和 `check_packages.R`：去重、英文日志、核心包与可选引擎包分离；可选包失败不会阻止 GUI 启动。
- 修复 `config.yml` 默认参数：随机效应默认改为 `none`，坐标列改为示例数据中的 `x`/`y`，最小示例不再强制要求 `studyDesign.csv`。
- 修复快速开始文档：示例 Y 文件改为真实存在的 `examples/Y.csv`，示例公式改为 `~ pH + moisture + canopy`。

## 模型运行状态

- HMSC：真实拟合已连接。默认 `hmsc_real_fit = TRUE` 时会执行 S1-S7：定义模型、MCMC 拟合、收敛诊断、模型拟合、模型拟合图、参数估计、环境梯度预测。需要 R 包 `Hmsc`，可选 `coda`、`corrplot`。
- jSDM：当前 GUI 分支为 scaffold/diagnostic 输出。项目内有 `R/engines/engine_jsdm.R` 真实适配器代码，但 GUI 尚未接入该适配器。
- GJAM：当前为 scaffold/diagnostic 输出。真实拟合需要接入 `gjam` production adapter。
- spOccupancy：当前为 scaffold/diagnostic 输出。真实拟合需要接入 `spOccupancy` production adapter。
- sjSDM：当前为 scaffold/diagnostic 输出。真实拟合需要 `sjSDM`、`reticulate`、Python/PyTorch，GPU 模式还需要 CUDA。
- boral：当前为 scaffold/diagnostic 输出。真实拟合需要 `boral`、`R2jags`，并且 Windows 另需安装 JAGS。

## 最小测试建议

- HMSC：上传 `examples/Y.csv` 与 `examples/XData.csv`，`distr = probit`，`XFormula = ~ .` 或 `~ pH + moisture + canopy`，`random effect = none`，`samples = 20`，`transient = 10`，`thin = 1`，`nChains = 2`，`nParallel = 1`。
- jSDM：上传 `examples/Y.csv` 与 `examples/XData.csv`，选择 `binomial_probit`；当前会生成 scaffold ZIP 与 diagnostics。
- GJAM：上传 `examples/Y.csv` 与 `examples/XData.csv`，`single typeName = PA`；当前会生成 scaffold ZIP 与 diagnostics。
- spOccupancy：需要检测/未检测结构数据；用普通 HMSC 示例只适合检查失败/诊断 ZIP 测试，不适合真实 occupancy 解释。
- sjSDM：上传 `examples/Y.csv` 与 `examples/XData.csv`，family 选 `binomial_probit`，device 选 `cpu`；当前会生成 scaffold ZIP 与 diagnostics。
- boral：上传 `examples/Y.csv` 与 `examples/XData.csv`，family 选 `binomial`；当前会生成 scaffold ZIP 与 diagnostics。

## 失败时先看哪里

- `diagnostics/engine_status.json`
- `diagnostics/data_check_messages.csv`
- `diagnostics/session_info.txt`
- `diagnostics/*_error.txt`
- `tables/data_summary.csv`
- `standard/run_summary.csv`

## 仍存在风险

- 只有 HMSC 已在 GUI 中真实执行拟合；其他模型仍需把 production adapter 接进 GUI run 函数。
- 当前环境未检测到 `Rscript.exe`，因此本次无法实际执行 R 语法检查或最小 HMSC 拟合测试。
- `external_packages/` 中包含第三方包源码、PDF、mhtml 和 zip，仅做资料/离线参考；生产安装仍建议使用 CRAN 或官方源。
- 旧 backup 文件仍保留，可能包含过时文字或历史逻辑；主入口以当前 `app.R`、`main.R`、`run_shiny_backend.R` 为准。
