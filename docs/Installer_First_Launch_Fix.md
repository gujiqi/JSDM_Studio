# 安装包首次启动修复说明

如果安装后的 JSDM Studio 弹出：

```text
Shiny backend did not start
R/Shiny backend did not start in time
```

常见原因：

1. 第一次启动正在安装 R 包，耗时超过等待时间；
2. R 包缺失；
3. 网络连接导致 R 包安装失败；
4. startup_log.txt 中有具体 R 报错。

本修复版做了：

1. `run_shiny_backend.R` 启动时自动检查并安装缺失 R 包；
2. 启动器等待时间从 5 分钟提高到 15 分钟；
3. 安装包开始菜单增加：
   - Repair R packages
   - Open startup log
   - Launch in default browser
4. 如果失败，可以直接打开 startup_log.txt 排查。

## 用户遇到问题时

让用户按顺序做：

```text
1. 开始菜单 → JSDM Studio → Repair R packages
2. 等待安装完成
3. 再打开 JSDM Studio
4. 如果仍失败，开始菜单 → Open startup log
5. 把 startup_log.txt 发给开发者
```
