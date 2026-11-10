# 构建与运行

[English](BUILD-AND-RUN.md) · [简体中文](BUILD-AND-RUN.zh-CN.md)

使用 Linux x86-64、GCC、支持 C++17 和 OpenMP 的 G++、Python 3.9 或更新版本、兼容的 NVIDIA GPU 驱动、CUDA 12.9 开发文件和 PanguLU 26.46 SDK。程序在单进程、单 GPU 上运行。SDK 和 GPU 许可证可向[解得科技](https://solved-tech.com/)获取。

在项目目录设置实际安装路径，按 SDK 说明配置许可证，再执行以下命令。

```bash
export PANGULU_ROOT=/path/to/pangulu-26.46
export CUDA_HOME=/usr/local/cuda-12.9

python3 tools/build.py --pangulu "$PANGULU_ROOT" --cuda "$CUDA_HOME"

for h in 85 75 65; do
    python3 tools/run_case.py "marmousi3d-h${h}" --output "_work/runs/marmousi3d-h${h}"
done
```

默认构建目录为 `_work/build`。可通过构建命令的 `--output` 和运行命令的 `--build` 指定其他目录。若 CUDA 开发文件分开安装，`--cuda` 指向运行库开发目录，`--cuda-nvcc-include` 指向附加头文件目录。构建器为各网格写出动态库检查日志。集群用户应在分配到 GPU 的作业内运行，并按 SDK 说明通过 `LD_LIBRARY_PATH` 提供所需运行库。

每次运行使用新的输出目录。运行器从真值模型生成 3 Hz 和 5 Hz 观测，再用平滑初始模型执行反演。使用 `--prepare-only` 可只复制案例配置和原始输入。

| 文件 | 内容 |
| --- | --- |
| `case.json` | 网格、频率、矩阵阶数和存储项数 |
| `observations/vp-true.bin`、`observations/vp-initial.bin`、`vp-final.bin` | 真值、初始和最终纵波速度。文件为小端 float64 节点值，单位 m/s，x 下标变化最快，无文件头。 |
| `observations/data-{0,1}.bin`、`predicted-{0,1}.bin` | 3 Hz 和 5 Hz 下接收器的三分量复位移，按接收器排序，每分量依次存储小端 float64 实部和虚部。 |
| `generate.*.log`、`invert.*.log`、`validation.json` | 生成与反演日志、验证结果、迭代轨迹和原生阶段计时 |

运行器核对每个原始 COO 系统的有限相对残差（≤ 1e-8）、正演接收器记录、有限差分梯度检查、12 次速度更新及波形和速度误差。核心耗时累计正式反演中的分析、数值分解和求解。合成观测、梯度检查、额外正演验证和应用墙钟分别记录或统计。
