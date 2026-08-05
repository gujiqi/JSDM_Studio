# 下一步：从 JSDM Studio 到安装包 / exe 启动器

第一版建议先用 `run_app.bat`。等 GUI 稳定后，再做安装包。

## 推荐路线

1. JSDM Studio + run_app.bat：当前版本，最稳。
2. Inno Setup 安装包：生成 `JSDM Studio_Setup.exe`，安装后桌面有图标。
3. exe 启动器：启动器只负责调用 Rscript 运行 Shiny app。
4. 便携 R / Electron：让用户基本感觉不到 R，但维护难度最高。

## 不建议第一步做的事情

不要尝试把 Hmsc R 包真正编译成单文件 exe。Hmsc 是 R 包，依赖 R 和多个 R 包；科研软件更适合做图形化工作流和安装包。
