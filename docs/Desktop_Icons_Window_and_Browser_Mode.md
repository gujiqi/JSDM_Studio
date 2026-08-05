# 桌面图标设置说明：窗口模式 + 浏览器模式

本版本安装时会默认勾选桌面图标，并创建两个入口。

## 桌面图标 1：JSDM Studio

默认主图标，打开独立窗口模式：

```text
JSDM Studio
→ JSDMStudioLauncher.exe
→ WebView2 独立窗口
```

适合普通用户，看起来更像 Windows 软件。

## 桌面图标 2：JSDM Studio Browser Mode

备用图标，打开浏览器模式：

```text
JSDM Studio Browser Mode
→ Launch_in_default_browser.bat
→ 默认浏览器打开 Shiny
```

适合独立窗口出问题时使用，稳定性最高。

## 开始菜单入口

安装后开始菜单会有：

```text
JSDM Studio - Window Mode
JSDM Studio - Browser Mode
Repair R packages
Open startup log
Troubleshoot R
Uninstall JSDM Studio
```

## 安装时默认勾选

Inno Setup 里已经去掉 `Flags: unchecked`，所以桌面图标默认勾选。

如果你不想给用户两个桌面图标，可以在安装时取消 Browser Mode 图标，或者后续改成只在开始菜单保留 Browser Mode。
