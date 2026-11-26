# 动态端口修复版

你日志里的核心错误是：

```text
createTcpServer: address already in use
Error in initialize(...) : Failed to create server
```

原因：固定端口 `3838` 已经被另一个 Shiny/R 进程占用。旧启动器检测到 3838 有服务，就误以为自己的后台启动成功；随后新的 R 后台再绑定 3838 就失败了。

本版本修复：

1. 启动器每次自动向 Windows 申请一个空闲端口；
2. 把这个端口传给 `run_shiny_backend.R`；
3. WebView2 也打开同一个动态端口；
4. 如果 R 进程提前退出，等待阶段会直接失败，不再误判为 ready。

## 运行

1. 双击 `Build_WebView2_Launcher.bat`
2. 双击 `Debug_WebView2_Launcher.bat`
3. 如果正常，再用 `Start_WebView2_Window.bat`

如果还有问题，把新的 `startup_log.txt` 发给我。
