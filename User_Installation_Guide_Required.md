# JSDM Studio 用户安装说明（普通用户版）

## 一句话说明

JSDM Studio 是一个基于 **R + Shiny + Hmsc** 的图形化联合物种分布模型（JSDM）工具。

用户看到的是一个 Windows 独立窗口，但后台仍然需要 R 和 Hmsc。

---

# 1. 用户只装我的 exe 就行吗？

## 不能完全只装你的 exe。

当前版本是：

```text
方案 B3：
JSDMStudio_Setup.exe
+ WebView2 独立窗口
+ 后台 R/Shiny/Hmsc
```

所以用户电脑还需要：

```text
1. R for Windows
2. Microsoft Edge WebView2 Runtime
3. Hmsc / shiny 等 R 包
```

但是：

```text
R 包不需要用户一个一个手动装。
你的程序会检查并尝试自动安装。
```

也就是说，普通用户实际操作可以写成：

```text
1. 先安装 R for Windows
2. 双击 JSDMStudio_Setup.exe
3. 双击桌面 JSDM Studio 图标
4. 第一次启动时程序自动检查/安装 R 包
```

如果电脑已有 WebView2 Runtime，则不需要额外安装 WebView2。

---

# 2. 用户必须安装什么？

## 必须安装 1：R for Windows

JSDM Studio 必须依赖 R，因为 Hmsc 是 R 包。

用户需要安装：

```text
R for Windows x64
```

官方下载：

```text
https://cran.r-project.org/bin/windows/base/
```

建议版本：

```text
R 4.3 或更高
推荐 R 4.4 / R 4.5
```

安装 R 时一路 Next 即可。

安装后，可以不用打开 RStudio。

---

## 必须安装 2：Microsoft Edge WebView2 Runtime

JSDM Studio 的独立窗口需要 WebView2 Runtime。

大多数 Windows 10 / Windows 11 已经自带。

如果打不开独立窗口，或者提示 WebView2 错误，用户需要安装：

```text
Microsoft Edge WebView2 Evergreen Runtime
```

官方下载：

```text
https://developer.microsoft.com/microsoft-edge/webview2/
```

选择：

```text
Evergreen Runtime
```

---

# 3. 用户不需要安装什么？

普通用户不需要安装：

```text
RStudio
Git
Node.js
Electron
Visual Studio
Inno Setup
.NET SDK
```

这些是开发者才需要的。

用户只需要：

```text
R
JSDMStudio_Setup.exe
```

以及在少数情况下安装 WebView2 Runtime。

---

# 4. JSDM Studio 需要哪些 R 包？

JSDM Studio 会检查以下 R 包：

```r
shiny
bslib
DT
yaml
ggplot2
zip
Hmsc
coda
ape
corrplot
htmltools
```

这些包的作用：

| R 包 | 作用 |
|---|---|
| `Hmsc` | 核心建模包，运行 Hmsc/JSDM |
| `shiny` | 图形界面后台 |
| `bslib` | 界面主题和布局 |
| `DT` | 表格预览 |
| `yaml` | 保存参数配置 used_config.yml |
| `ggplot2` | 出图 |
| `zip` | 打包结果 |
| `coda` | MCMC 诊断 |
| `ape` | 系统发育树 |
| `corrplot` | 物种关联图 |
| `htmltools` | HTML 报告 |

---

# 5. R 包怎么安装？

## 推荐方式：让 JSDM Studio 自动安装

用户第一次启动时，程序会运行：

```text
check_packages.R
```

它会检查哪些包缺失。

如果缺包，它会调用：

```text
install_packages.R
```

自动安装。

用户只需要保持网络畅通。

---

## 手动安装方式

如果自动安装失败，用户可以进入 JSDM Studio 安装目录，双击：

```text
install_packages.bat
```

也可以打开 R，运行：

```r
options(repos = c(CRAN = "https://cloud.r-project.org"))

packages <- c(
  "shiny", "bslib", "DT", "yaml", "ggplot2", "zip",
  "Hmsc", "coda", "ape", "corrplot", "htmltools"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}
```

---

# 6. 用户安装顺序

## 最推荐顺序

```text
第 1 步：安装 R for Windows
第 2 步：双击 JSDMStudio_Setup.exe
第 3 步：桌面出现 JSDM Studio 图标
第 4 步：双击 JSDM Studio
第 5 步：第一次启动时等待 R 包检查/安装
第 6 步：独立窗口打开
第 7 步：上传数据并运行模型
```

---

# 7. 如果打不开怎么办？

## 情况 1：提示找不到 R

说明用户没有安装 R，或者 R 没有安装在常见位置。

解决：

```text
安装 R for Windows
https://cran.r-project.org/bin/windows/base/
```

安装后重新打开 JSDM Studio。

---

## 情况 2：提示缺 R 包

进入 JSDM Studio 文件夹，双击：

```text
install_packages.bat
```

或者在开始菜单里找：

```text
Repair packages
```

如果有这个入口。

---

## 情况 3：WebView2 窗口打不开

安装：

```text
Microsoft Edge WebView2 Evergreen Runtime
```

下载：

```text
https://developer.microsoft.com/microsoft-edge/webview2/
```

---

## 情况 4：独立窗口有问题

使用备用浏览器模式：

```text
Launch_in_default_browser.bat
```

这个模式会用普通浏览器打开界面，稳定性最高。

---

## 情况 5：启动后报错

查看日志文件：

```text
startup_log.txt
```

这个文件在 JSDM Studio 安装目录中。

把最后 30 行发给开发者。

---

# 8. 推荐给用户的最短说明

你可以直接发给用户这段：

```text
使用前请先安装 R for Windows：
https://cran.r-project.org/bin/windows/base/

然后双击 JSDMStudio_Setup.exe 安装。
安装后桌面会出现 JSDM Studio 图标。

第一次启动时，程序会自动检查并安装所需 R 包：
Hmsc、shiny、bslib、DT、yaml、ggplot2、zip、coda、ape、corrplot、htmltools。

大多数 Windows 10/11 已经自带 WebView2 Runtime。
如果独立窗口打不开，请安装 Microsoft Edge WebView2 Evergreen Runtime：
https://developer.microsoft.com/microsoft-edge/webview2/
```

---

# 9. 当前版本定位

当前版本是：

```text
方案 B3：安装包版 + WebView2 独立窗口
```

它不是完全不依赖 R 的商业软件。

它的优势是：

```text
用户不用写 R 代码
用户不用打开 RStudio
用户不用手动找 app.R
用户通过桌面图标打开独立窗口
```

它仍然依赖：

```text
R
Hmsc R 包
Shiny
WebView2 Runtime
```

未来如果要做到“用户只装一个 exe，连 R 都不用装”，需要做：

```text
portable R
+ 预装 R 包库
+ WebView2 Runtime 检测/安装
+ 更大的安装包
```
