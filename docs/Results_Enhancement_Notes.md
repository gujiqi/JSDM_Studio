# 结果增强说明

本版本根据 S1-S7 的脚本式工作流扩展输出：

- S1：定义未拟合模型，保存 unfitted_models.RData
- S2：拟合模型，保存 models_thin_*_samples_*_chains_*.RData
- S3：收敛诊断，保存 MCMC_convergence.txt 和 MCMC_traceplots.pdf
- S4/S5：模型拟合与交叉验证，保存 model_fit 文件和图
- S6：参数估计，保存 Beta/Gamma/Omega/方差分解表和图
- S7：环境梯度预测，保存 predictions_environmental_gradients.pdf

注意：GUI 是单模型版本，不是官方脚本里多模型批处理的完整复制。后续可以继续做“多模型比较版”。
