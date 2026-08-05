# 开发者打包说明：方案 B3 WebView2 版

## 需要安装

开发者电脑需要：

```text
R for Windows
Inno Setup
.NET SDK
```

如果使用 .NET 5 兼容版，可以用：

```text
.NET 5 SDK
```

如果使用更推荐的新版本，可以用：

```text
.NET 8 SDK x64
```

普通用户不需要 .NET SDK。

---

## 开发者构建顺序

进入：

```text
JSDMStudio
```

先构建 WebView2 启动器：

```text
Build_WebView2_Launcher.bat
```

测试独立窗口：

```text
Debug_WebView2_Launcher.bat
```

确认成功后，生成 Setup.exe：

```text
JSDMStudio\installer\build_installer.bat
```

输出位置：

```text
JSDMStudio\installer\output\JSDMStudio_Setup.exe
```

---

## 发布给用户

只发：

```text
JSDMStudio_Setup.exe
```

同时告诉用户：

```text
需要先安装 R for Windows
```

如果窗口打不开，再安装：

```text
Microsoft Edge WebView2 Runtime
```

---

## 当前架构

```text
JSDMStudioLauncher.exe
→ 找 Rscript.exe
→ 启动 run_shiny_backend.R
→ Shiny 在本地端口运行
→ WebView2 窗口加载本地界面
```

这不是 Electron，不需要 Node.js。
