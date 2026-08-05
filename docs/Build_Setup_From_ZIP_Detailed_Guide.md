# 从这个 ZIP 生成 JSDMStudio_Setup.exe：完整步骤

这个说明给开发者/作者使用。普通用户不需要看这个文件。

你现在拿到的是：

```text
JSDMStudio_WebView2_B3_NET5_HIGHDPI_WITH_DOCS.zip
```

它还不是最终发给用户的安装包。  
你需要先在自己电脑上把它编译成：

```text
JSDMStudio_Setup.exe
```

然后把这个 `JSDMStudio_Setup.exe` 发给普通用户。

---

# 一、你需要先安装什么？

开发者电脑需要安装：

## 1. R for Windows

下载：

```text
https://cran.r-project.org/bin/windows/base/
```

你已经有 R 就不用重复装。

---

## 2. Inno Setup

用来生成 `JSDMStudio_Setup.exe` 安装包。

下载：

```text
https://jrsoftware.org/isinfo.php
```

安装后，电脑里应该有：

```text
ISCC.exe
```

---

## 3. .NET SDK

用来编译 WebView2 独立窗口启动器：

```text
JSDMStudioLauncher.exe
```

你现在如果是 .NET 5 SDK，可以用本 ZIP 的 NET5 兼容版。

检查：

```powershell
dotnet --version
dotnet --list-sdks
```

如果显示：

```text
5.0.408
```

也可以继续。

如果你想用新版本，推荐安装：

```text
.NET 8 SDK x64
```

下载：

```text
https://dotnet.microsoft.com/download
```

---

# 二、生成安装包的最短流程

解压 ZIP 后，你会看到：

```text
JSDMStudio/
```

按这个顺序运行：

```text
1. 进入 JSDMStudio
2. 双击 Build_WebView2_Launcher.bat
3. 双击 Start_WebView2_Window.bat 测试窗口
4. 如果测试成功，进入 JSDMStudio\installer
5. 双击 build_installer.bat
6. 去 JSDMStudio\installer\output 找 JSDMStudio_Setup.exe
```

一句话：

```text
Build_WebView2_Launcher.bat
→ Start_WebView2_Window.bat
→ installer\build_installer.bat
→ installer\output\JSDMStudio_Setup.exe
```

---

# 三、详细步骤

## 第 1 步：解压 ZIP

把 ZIP 解压到一个简单路径，建议不要放在太深或带特殊符号的位置。

推荐：

```text
C:\Users\你的用户名\Downloads\JSDMStudio_Build\
```

不要放在：

```text
OneDrive 同步目录
中文特别复杂的路径
很深的多层目录
```

---

## 第 2 步：进入 JSDMStudio 文件夹

打开：

```text
JSDMStudio/
```

你应该看到：

```text
app.R
R/
Build_WebView2_Launcher.bat
Start_WebView2_Window.bat
Debug_WebView2_Launcher.bat
install_packages.bat
installer/
webview2_launcher/
```

---

## 第 3 步：先安装 R 包

双击：

```text
install_packages.bat
```

它会安装/检查：

```text
Hmsc
shiny
bslib
DT
yaml
ggplot2
zip
coda
ape
corrplot
htmltools
```

如果已经装过，会显示已经存在。

---

## 第 4 步：编译 WebView2 启动器

双击：

```text
Build_WebView2_Launcher.bat
```

成功后，`JSDMStudio` 文件夹里应该出现：

```text
JSDMStudioLauncher.exe
JSDMStudioLauncher.dll
JSDMStudioLauncher.deps.json
JSDMStudioLauncher.runtimeconfig.json
Microsoft.Web.WebView2.Core.dll
Microsoft.Web.WebView2.WinForms.dll
Microsoft.Web.WebView2.Wpf.dll
```

这些文件都要保留。

如果只看到 `JSDMStudioLauncher.exe`，没有 DLL/json，说明依赖没有复制完整。

---

## 第 5 步：测试独立窗口

双击：

```text
Start_WebView2_Window.bat
```

如果它打开一个独立 Windows 窗口，并显示 JSDM Studio 界面，说明窗口启动器成功。

如果有问题，双击：

```text
Debug_WebView2_Launcher.bat
```

然后查看：

```text
startup_log.txt
```

---

## 第 6 步：生成 Setup.exe

确认窗口能打开后，进入：

```text
JSDMStudio\installer
```

双击：

```text
build_installer.bat
```

它会调用 Inno Setup 生成安装包。

成功时会看到：

```text
Successful compile
Build completed
```

---

## 第 7 步：找到最终安装包

安装包位置：

```text
JSDMStudio\installer\output\
```

文件名：

```text
JSDMStudio_Setup.exe
```

这个文件就是最终可以发给用户的安装包。

---

# 四、你发给用户什么？

普通用户只需要这个：

```text
JSDMStudio_Setup.exe
```

不要把整个开发 ZIP 发给普通用户。

同时告诉用户：

```text
请先安装 R for Windows：
https://cran.r-project.org/bin/windows/base/
```

如果独立窗口打不开，再安装：

```text
Microsoft Edge WebView2 Evergreen Runtime：
https://developer.microsoft.com/microsoft-edge/webview2/
```

---

# 五、用户安装后会看到什么？

用户双击：

```text
JSDMStudio_Setup.exe
```

安装后桌面出现：

```text
JSDM Studio
```

用户双击桌面图标后：

```text
JSDMStudioLauncher.exe
→ 自动寻找 Rscript.exe
→ 启动 run_shiny_backend.R
→ WebView2 独立窗口打开 JSDM Studio
```

用户不需要打开 RStudio，也不需要写 R 代码。

---

# 六、如果 build_installer.bat 报错怎么办？

## 报错：找不到 Inno Setup / ISCC.exe

说明你没装 Inno Setup，或者 PATH 没有识别。

安装：

```text
https://jrsoftware.org/isinfo.php
```

然后重新双击：

```text
installer\build_installer.bat
```

---

## 报错：找不到 dotnet

说明没有安装 .NET SDK。

检查：

```powershell
dotnet --version
```

安装 .NET SDK 后重试。

---

## 报错：WebView2 启动器编译失败

先单独运行：

```text
Build_WebView2_Launcher.bat
```

看具体错误。

如果是 .NET 版本不匹配，可以使用本包的 NET5 兼容版本，或安装 .NET 8 SDK x64。

---

## 报错：R/Shiny 后台启动失败

看：

```text
startup_log.txt
```

常见原因：

```text
R 包没装全
端口被占用
R 路径不对
app.R 有错误
```

---

# 七、完整发布建议

正式发布时建议准备三个文件：

```text
1. JSDMStudio_Setup.exe
2. 用户安装最短说明.txt
3. 示例数据说明.pdf 或 README
```

给用户的信息可以写：

```text
请先安装 R for Windows，然后运行 JSDMStudio_Setup.exe。
第一次启动会自动检查/安装 Hmsc、shiny 等 R 包。
如果独立窗口打不开，请安装 Microsoft Edge WebView2 Runtime 或使用备用浏览器模式。
```

---

# 八、当前版本定位

当前版本是：

```text
方案 B3：
Inno Setup 安装包
+ WebView2 独立窗口启动器
+ R/Shiny/Hmsc 后台
```

它不是完全独立、不依赖 R 的商业软件。

它的优点是：

```text
用户不用写代码
用户不用找 app.R
用户不用打开 RStudio
用户通过桌面图标打开独立窗口
```

它仍然需要：

```text
R for Windows
R 包
WebView2 Runtime
```

未来如果要“用户只装一个 exe，连 R 都不用装”，需要继续做：

```text
portable R
+ 预装 R 包库
+ WebView2 Runtime 检测/安装
+ 更大的安装包
```
