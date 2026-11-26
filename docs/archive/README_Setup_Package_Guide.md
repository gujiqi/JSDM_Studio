# JSDM Studio Setup.exe 安装包版说明

这个文件夹是方案 B：Windows 安装包版。

## 用户看到什么

用户最终会得到：

- `JSDMStudio_Setup.exe`
- 桌面图标：`JSDM Studio`
- 开始菜单：`JSDM Studio`

用户双击桌面图标后，程序会：

1. 自动寻找 Rscript.exe；
2. 如果没有 R，会打开 R 官方下载页；
3. 检查 Hmsc、shiny、yaml 等 R 包；
4. 缺包时自动运行 `install_packages.R`；
5. 启动 Shiny 图形界面。

## 你如何生成 Setup.exe

你的电脑需要先安装 Inno Setup：

https://jrsoftware.org/isinfo.php

安装后，进入本文件夹：

```text
JSDMStudio/installer/
```

双击：

```text
build_installer.bat
```

成功后会生成：

```text
JSDMStudio_Setup.exe
```

位置通常在：

```text
JSDMStudio/installer/output/
```

## 注意

这个 Setup.exe 不是把 Hmsc 真正编译成机器码。它是安装和启动器：

```text
Setup.exe
→ 安装 JSDM Studio 文件
→ 创建桌面图标
→ 图标调用 Launch_JSDMStudio.bat
→ bat 调用 Rscript
→ Rscript 启动 Shiny/Hmsc
```

所以用户电脑仍然需要 R。区别是：用户不需要写 R 代码，不需要找 app.R，也不需要手动从命令行启动。
