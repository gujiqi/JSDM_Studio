# JSDM Studio 美化注释版说明

这个版本主要改进了界面：

1. 增加顶部导航：开始、数据上传、参数设置、运行与结果、参数解释。
2. 增加每个关键参数的中文注释。
3. 增加数据概览卡片：Y 矩阵大小、XData 大小、物种数、运行状态。
4. 增加更清楚的运行区和输出文件说明。
5. 保留原来的后端逻辑：仍然调用 Hmsc、sampleMcmc、computePredictedValues、evaluateModelFit 等函数。

使用方式：
1. 双击 install_packages.bat 安装依赖。
2. 双击 run_app.bat 启动界面。
3. 第一次建议用示例数据：
   - Y：Y_presence_absence.csv
   - XData：XData.csv
   - distr：probit
   - samples：20
   - transient：10
   - thin：1
   - nChains：2
   - nParallel：1

注意：
- 这个版本还是方案 A：用户需要安装 R，但不需要写 R 代码。
- 正式论文分析不能用太小的 MCMC 参数，需要提高 samples/transient 并检查诊断。
