# 构建与运行

[English](BUILD-AND-RUN.md) · [总览](../README.zh-CN.md) · [接入说明](INTEGRATION.zh-CN.md)

## 准备环境

使用 Linux x86-64、支持 C++17 的编译器、CUDA GPU 和匹配设备的 PanguLU 26.46 许可证。CUDA 工具链需与 SDK 相匹配。原案例采用 CUDA 12.9 版 PanguLU 和 A100 PCIe 40 GB，百万阶作业申请了 100 GiB 主机内存，可作为资源配置参考。

使用 Python 3.11 和 NumPy 后端。以下命令均从解压后的 `FEALPy-with-Solved-Tech` 目录执行，将 FEALPy 上游源码和 Python 环境放在同级目录。无图形界面的 Linux 服务器也需安装 Gmsh Python 包依赖的系统 OpenGL/GLU 库。

```bash
python3.11 -m venv ../fealpy-venv
source ../fealpy-venv/bin/activate
git clone https://github.com/weihuayi/fealpy.git ../fealpy-upstream
git -C ../fealpy-upstream checkout --detach e033533bd8e24f796f9c5284685b20a98b4d11cb
python -m pip install numpy==1.26.4 scipy==1.16.2
env -u WITH_PANGULU python -m pip install -e ../fealpy-upstream
```

## 构建适配库

将示例安装路径替换为本机目录，按 SDK 的授权说明配置许可证。若文件保存在自定义位置，将环境变量 `PANGULU_LICENSE_PATHS` 设置为该文件的绝对路径。多卡环境可通过 `CUDA_VISIBLE_DEVICES` 选择 GPU，适配器打开可见设备 0。

```bash
export PANGULU_ROOT=/path/to/pangulu-26.46
export CUDA_HOME=/usr/local/cuda-12.9
export LD_LIBRARY_PATH="$PANGULU_ROOT/lib:$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
bash scripts/build.sh
```

构建生成 `build/libfealpy_pangulu_real.so` 和 `build/libfealpy_pangulu_complex.so`。若 CUDA 组件分开安装，可在构建前分别设置 `CUDA_INCLUDE`、`CUDA_NVCC_INCLUDE`、`CUDA_LIB`，指向运行时头文件、编译器头文件和运行时库目录。运行时将实际的库目录和 SDK 所需的依赖库加入动态库搜索路径。

## 运行与校验

以下命令使用小网格，各调用一次 GPU 求解器。去掉 `--small` 即采用案例说明中的百万阶网格。

```bash
python cases/elastic3d/run.py --small
python cases/maxwell2d/run.py --small
python cases/advdiff2d/run.py --small
```

保存三个默认网格的解场系数：

```bash
python cases/elastic3d/run.py --output output/elastic3d.npy
python cases/maxwell2d/run.py --output output/maxwell2d.npy
python cases/advdiff2d/run.py --output output/advdiff2d.npy
```

终端依次报告矩阵阶数、CSR 存储项数、原方程相对残差和案例对应的物理场误差。通过后打印 `validation=passed`。需要保留结果时，为每次运行指定独立的 `.npy` 文件；其中是一维 NumPy 数组，按原有限元自由度顺序排列。系数的物理含义分别见[线弹性](../cases/elastic3d/README.zh-CN.md)、[频域电磁场](../cases/maxwell2d/README.zh-CN.md)和[标量输运](../cases/advdiff2d/README.zh-CN.md)说明。
