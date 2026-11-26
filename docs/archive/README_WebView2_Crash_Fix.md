# WebView2 闪退修复说明

你遇到“生成了但闪退”，最可能原因是：

上一版 `Build_WebView2_Launcher.bat` 只复制了：

```text
JSDMStudioLauncher.exe
```

但 .NET / WebView2 程序还需要同目录下的依赖文件，例如：

```text
JSDMStudioLauncher.dll
JSDMStudioLauncher.deps.json
JSDMStudioLauncher.runtimeconfig.json
Microsoft.Web.WebView2.Core.dll
Microsoft.Web.WebView2.WinForms.dll
```

只复制 exe 会导致双击后马上关闭。

本版本修复：

```text
copy /Y launcher_build\*.* .\
```

也就是把发布目录中的所有依赖一起复制到 JSDMStudio 文件夹。

## 运行顺序

1. 解压本 ZIP
2. 进入 `JSDMStudio`
3. 双击 `Build_WebView2_Launcher.bat`
4. 确认 JSDMStudio 文件夹里有：
   - JSDMStudioLauncher.exe
   - JSDMStudioLauncher.dll
   - Microsoft.Web.WebView2.WinForms.dll
5. 双击 `Start_WebView2_Window.bat`

如果还闪退，双击：

```text
Debug_WebView2_Launcher.bat
```

它会尝试打开 `startup_log.txt`，方便定位原因。
