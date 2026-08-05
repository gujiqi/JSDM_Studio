# .NET 5 兼容修复版

你电脑当前显示：

```text
dotnet --version
5.0.408
```

上一版 WebView2 启动器目标是 `net8.0-windows`，所以 .NET 5 SDK 无法编译，报：

```text
NETSDK1045 当前 .NET SDK 不支持将 .NET 8.0 设置为目标
```

本版本已经改成：

```text
<TargetFramework>net5.0-windows</TargetFramework>
```

并把 `ApplicationConfiguration.Initialize()` 改成 .NET 5 可用写法。

所以你现在可以直接用现有 .NET 5.0.408 测试。

## 运行顺序

1. 解压本 ZIP
2. 进入 `JSDMStudio`
3. 双击 `Build_WebView2_Launcher.bat`
4. 编译成功后双击 `Start_WebView2_Window.bat`
5. 如果窗口能打开，再进入 `JSDMStudio\installer`
6. 双击 `build_installer.bat`
7. 安装包在 `installer\output\JSDMStudio_Setup.exe`

## 更推荐

长期更推荐安装 .NET 8 SDK，但为了你现在能跑，这个包兼容 .NET 5。
