# 构建与运行

[English](BUILD-AND-RUN.md)

准备 Linux x86-64、GNU C/C++/Fortran 编译器、Make、Expat 开发库、BLAS/LAPACK 和 Python 3.12 或更新版本。PanguLU 使用 26.44、ABI 2 和 CUDA 12.9 版交付包，需要与目标 NVIDIA GPU 对应的许可证。运行采用单 PHG 进程、单 GPU。

## 取得上游源码并构建

在本目录的项目根目录运行。将示例路径换成自己的安装位置。

```bash
python3 tools/prepare.py
export PANGULU_ROOT=/path/to/pangulu-26.44
export CUDA_HOME=/usr/local/cuda-12.9
python3 tools/build.py --pangulu "$PANGULU_ROOT" --cuda "$CUDA_HOME"
```

`prepare.py` 从官方站点下载 `upstream.lock.json` 中的两个归档，解包至 `_work/clean`。离线时可把同名官方归档放入 `_work/downloads`。如需另建干净目录，给准备、构建和运行命令传同一个 `--workspace`。构建日志在工作目录的 `logs/`，可执行文件在 `bin/`。

非标准 BLAS/LAPACK 路径可用 `--blas` 和 `--lapack` 传入链接参数，例如 `--blas="-L/path/to/openblas/lib -lopenblas"`。CUDA 拆分安装时，可用 `--cuda-include`、`--cuda-nvcc-include` 和 `--cuda-lib` 分别指定头文件和运行库目录。PanguLU 的其他共享库依赖按 SDK 的安装说明配置。

## 运行三个案例

在有 GPU 的节点上，从项目根目录执行。每次运行指定新的输出目录。

```bash
export PANGULU_LICENSE_PATHS=/path/to/your-gpu.license.key
python3 tools/run_case.py solid-adaptive-stress --output _work/runs/adaptive --vtk
python3 tools/run_case.py solid-static-p2 --output _work/runs/static-p2 --vtk
python3 tools/run_case.py interface-sphere-dg1 --output _work/runs/interface
```

默认运行完整案例，PanguLU 主机线程数默认为 40，可用 `--solver-threads` 调整。在已分配的 GPU 节点上可用 `--profile connection-test` 先跑小网格。固体位移场的 VTK 输出由 `--vtk` 启用，界面案例使用上游自带的 VTK 输出。

`run-report.json` 记录组装尺寸、残差、原生阶段和场文件。运行验收检查矩阵尺寸、求解状态、原 CSR 相对残差和结果场。详细说明见[接口与输出](INTEGRATION.zh-CN.md)。
