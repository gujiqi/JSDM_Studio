# Inno Setup 架构字段修复说明

你遇到的错误：

```text
Architecture identifier "x64compatiblecompatible" is invalid
```

原因是之前脚本多次自动替换，把：

```text
x64
```

先改成：

```text
x64compatible
```

又重复改成了：

```text
x64compatiblecompatible
```

本版本已修复为：

```text
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
```

运行：

```text
JSDMStudio\installer\build_installer.bat
```

即可重新生成安装包。
