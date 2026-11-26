# High-DPI 清晰修复版

这个版本针对 WebView2 窗口发虚/不清晰做了修复：

1. 程序启动时启用 PerMonitorV2 DPI awareness；
2. 增加 app.manifest 高 DPI 声明；
3. WinForms 窗口设置 `AutoScaleMode = Dpi`；
4. WebView2 设置 `ZoomFactor = 1.0`；
5. 默认字体改为 `Microsoft YaHei UI`，中文更清楚。

## 使用

1. 解压本 ZIP
2. 进入 `JSDMStudio`
3. 双击 `Build_WebView2_Launcher.bat`
4. 双击 `Debug_WebView2_Launcher.bat` 或 `Start_WebView2_Window.bat`

## 如果还不清楚

右键 `JSDMStudioLauncher.exe`：

```text
属性
→ 兼容性
→ 更改高 DPI 设置
→ 勾选“替代高 DPI 缩放行为”
→ 缩放执行：应用程序
```

但正常情况下，本版本已经不需要手动设置。
