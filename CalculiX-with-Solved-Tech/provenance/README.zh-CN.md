# 上游源码与案例来源

[English](README.md) | [简体中文](README.zh-CN.md) · [首页](../README.zh-CN.md)

`upstream/ccx_2.23.src.tar.bz2`包含完整的CalculiX2.23源码发行包，构建工具从该包解压并编译。`upstream/src/`中的部分源码便于查看求解接口及可选计时补丁。

来源：[CalculiX下载页面](https://www.dhondt.de/)。源码中保留版权、作者和GPL声明。

案例输入及include文件位于`cases/<case>/source/`。各例`parameters.json`记录来源URL、版本、网格参数及实测方程组规模。

PanguLU和CUDA作为外部依赖，其SDK与许可证单独安装配置。
