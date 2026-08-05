# C# 编译错误修复说明

你遇到的错误：

```text
error CS1514: 应为 {
error CS1513: 应输入 }
```

原因是上一版 `Program.cs` 使用了较新的 C# 文件范围命名空间写法：

```csharp
namespace JSDMStudioLauncher;
```

但你当前 .NET 5 SDK / MSBuild 环境不支持这个写法。

本版本已改成 .NET 5 兼容写法：

```csharp
namespace JSDMStudioLauncher
{
    ...
}
```

并把 `LangVersion` 固定为 `9.0`。

运行顺序：

1. 解压本 ZIP
2. 进入 `JSDMStudio`
3. 双击 `Build_WebView2_Launcher.bat`
4. 成功后双击 `Start_WebView2_Window.bat`
