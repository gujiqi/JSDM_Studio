# 点击按钮闪退调试版

这个版本加强了日志和错误捕获：

- 捕获 C# UI 异常
- 捕获 WebView2 渲染进程失败
- 监控 R/Shiny 后台是否退出
- 所有关键错误写入 `startup_log.txt`
- 增加 `Debug_WebView2_Launcher.bat`

## 运行

1. 双击 `Build_WebView2_Launcher.bat`
2. 双击 `Debug_WebView2_Launcher.bat`
3. 在窗口里点击你之前导致闪退的按钮
4. 如果还闪退，把 `startup_log.txt` 发给我

## 备用

如果 WebView2 版不稳定，先用：

```text
Launch_in_default_browser.bat
```

浏览器版稳定性最高，适合作为正式发布 fallback。
