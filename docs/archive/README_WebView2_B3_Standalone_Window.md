# JSDM Studio 方案 B3：WebView2 独立窗口启动器

这是你要的折中方案：

```text
方案 B 安装包
+ 独立 Windows 窗口
+ 不打开 Edge/Chrome/默认浏览器
+ 后台仍然用 R/Shiny/Hmsc
```

## 它怎么工作？

```text
JSDMStudioLauncher.exe
→ 自动寻找 Rscript.exe
→ 启动 run_shiny_backend.R
→ 后台运行 Shiny: http://127.0.0.1:3838
→ WebView2 窗口加载这个本地界面
```

用户看到的是 Windows 窗口，不是浏览器标签页。

## 用户需要额外安装什么？

用户电脑需要：

1. R for Windows
2. Hmsc/Shiny 等 R 包
3. Microsoft Edge WebView2 Runtime

Windows 10/11 大多数电脑已经有 WebView2 Runtime。  
如果没有，程序会提示安装。

## 开发者需要什么？

要编译 `JSDMStudioLauncher.exe`，你的电脑需要：

```text
.NET 8 SDK
```

下载地址：

```text
https://dotnet.microsoft.com/download
```

## 测试 WebView2 窗口

进入 JSDMStudio 文件夹，双击：

```text
Start_WebView2_Window.bat
```

第一次如果没有 launcher，会自动编译。

## 生成 Setup.exe

进入：

```text
JSDMStudio\installer
```

双击：

```text
build_installer.bat
```

它会先编译 WebView2 启动器，再用 Inno Setup 生成安装包。

输出：

```text
JSDMStudio\installer\output\JSDMStudio_Setup.exe
```

## 备用入口

如果 WebView2 出问题，开始菜单里保留：

```text
Launch in default browser
Troubleshoot R
```

这样用户不会完全卡死。
