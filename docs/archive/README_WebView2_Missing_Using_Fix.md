# WebView2 启动器 using 修复版

你遇到的错误：

```text
CS0246: 未能找到类型或命名空间名 Form / EventArgs / FormClosedEventArgs / Version / TimeSpan / STAThread
```

原因是 .NET 5 / C# 9 环境没有自动启用较新的 implicit usings，所以 Program.cs 里缺少显式 using。

本版本已在 Program.cs 顶部加入：

```csharp
using System;
using System.Windows.Forms;
using System.Drawing;
using System.Diagnostics;
using System.Net.Sockets;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading;
using Microsoft.Web.WebView2.WinForms;
```

并把 csproj 里的 implicit usings 关闭，固定为 C# 9。

运行顺序：

1. 解压本 ZIP
2. 进入 `JSDMStudio`
3. 双击 `Build_WebView2_Launcher.bat`
4. 成功后双击 `Start_WebView2_Window.bat`
5. 如果窗口成功，再去 `installer` 里运行 `build_installer.bat`
