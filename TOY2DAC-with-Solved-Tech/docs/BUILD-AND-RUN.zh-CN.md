# 构建与运行

[English](BUILD-AND-RUN.md) · [简体中文](BUILD-AND-RUN.zh-CN.md)

使用 Linux x86-64、GCC、支持 C++17 的 G++、Python 3.9 或更新版本、兼容的 NVIDIA GPU 驱动、CUDA 12.9 开发文件和 PanguLU 26.46 SDK。三维程序以单进程运行。SDK 和适用于目标 GPU 的许可证可向[解得科技](https://solved-tech.com/)获取。

在项目根目录运行下面的命令，并替换为实际安装路径。按 SDK 的说明配置许可证。本次性能测试使用的默认授权文件位置为 `~/.pangulu/license.key`。

```bash
export PANGULU_ROOT=/path/to/pangulu-26.46
export CUDA_HOME=/usr/local/cuda-12.9

python3 tools/build.py --pangulu "$PANGULU_ROOT" --cuda "$CUDA_HOME"

for h in 85 75 65; do
    python3 tools/run_case.py "marmousi3d-h${h}" --output "_work/runs/marmousi3d-h${h}"
done
```

默认构建目录是 `_work/build`。如需改用其他目录，构建时传 `--output`，运行时传 `--build`。若 CUDA 开发文件分开安装，`--cuda` 指向运行库开发目录，`--cuda-nvcc-include` 指向附加头文件目录。构建和运行时，通过 `LD_LIBRARY_PATH` 提供 SDK 所需的 CUDA 运行库。构建器将动态库检查写入 `loader-marmousi3d_pangulu.txt`。在集群上请在已分配的 GPU 作业内运行。

每次运行使用新的输出目录。加 `--prepare-only` 可单独准备模型数据和案例参数，无需先构建或使用 GPU。

| 文件 | 内容 |
| --- | --- |
| `case.json`、`metadata.txt` | 案例参数、网格尺寸和应用计算的残差 |
| `model.vtk` | 三维纵波和横波速度，单位 m/s |
| `field.vtk` | 位移向量的实部、虚部和合成幅值 |
| `field_complex.bin` | 完整三维位移。开头为三个小端 32 位节点数，其后每个节点依次写入 x、y、z 三分量的实部和虚部，共六个 float64。x 下标变化最快。 |
| `midplane_y.csv` | 横向中间剖面上的速度和复位移 |
| `stdout.log`、`stderr.log`、`validation.json` | 应用与求解器日志、残差、波场检查及原生阶段计时 |

运行入口核对矩阵阶数和 COO 存储项数，检查原始方程的相对残差有限且不超过 `1e-8`，并验证完整二进制位移的尺寸、长度和数值有限性。核心耗时由原生分析、数值分解和求解相加，建模、组装、显式数据传输、输出和校验在核心计时之外。发生求解错误时，可从日志查看出错阶段和状态码。
