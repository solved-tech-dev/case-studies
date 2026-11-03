# 构建与运行

[English](BUILD-AND-RUN.md) | [简体中文](BUILD-AND-RUN.zh-CN.md) · [首页](../README.zh-CN.md)

以下命令均在`Calculix/`目录中运行。使用集群时，请先申请GPU计算资源。

## 1. 环境依赖

| 组件 | 要求 |
|---|---|
| 平台 | Linux x86_64 |
| 构建工具 | GCC、G++、GFortran、GNU Make、ar、patch；Python 3.12或更新版本 |
| PanguLU | 26.42、ABI2、CUDA12.9，包含SDK头文件和动态库 |
| GPU运行环境 | 匹配的NVIDIA驱动及CUDA runtime开发文件 |
| 数值库 | LP64 BLAS/LAPACK、顺序ARPACK、SPOOLES2.2头文件及顺序库 |
| 其他库 | zlib开发文件 |
| 网格生成 | NumPy，仅`tools/generate_mesh.py`需要 |

参考运行配置为一块A10040GB GPU、40个主机CPU核和100GiB主机内存。运行工具要求至少20GiB可用磁盘空间。

CalculiX的约束处理使用SPOOLES辅助例程，相关模块的构建依赖ARPACK。主线性方程求解入口由PanguLU适配器提供。

如需购买PanguLU并获取许可证，请访问[解得科技官网](https://solved-tech.com/)。请单独安装SDK与许可证，并核对上述接口版本。

## 2. 编译CalculiX与适配器

将示例路径替换为本机依赖目录：

```bash
python3 tools/build.py \
  --pangulu-root /path/to/external/pangulu \
  --cuda-root /path/to/external/cuda \
  --spooles-root /path/to/external/SPOOLES.2.2 \
  --arpack-libs "/path/to/external/libarpack.a" \
  --blas-libs "-L/path/to/external/openblas/lib -lopenblas" \
  --workspace _work/build --jobs 16 --timing
```

脚本将`provenance/upstream/ccx_2.23.src.tar.bz2`解压到工作目录，编译CalculiX与适配器，生成`_work/build/ccx_pangulu`。构建采用默认32位CalculiX整数接口和LP64数值库。

| 选项 | 用途 |
|---|---|
| `--pangulu-root` | 包含`include/pangulu.h`及`lib/`下PanguLU库的目录 |
| `--cuda-root` | 包含`include/`及`lib64/`或`lib/`的CUDA目录 |
| `--spooles-root` | 包含SPOOLES头文件与`spooles.a`的目录 |
| `--arpack-libs` | 顺序ARPACK的链接选项或静态库路径 |
| `--blas-libs` | LP64 BLAS/LAPACK链接选项 |
| `--jobs` | 并行编译进程数 |
| `--timing` | 输出Newton各阶段计时 |
| `--resume` | 使用相同依赖路径和选项继续构建 |

使用拆分安装的CUDA时，可重复传入`--cuda-include`添加头文件目录，并用`--cuda-lib`指定runtime库目录。使用oneMKL时，通过`--blas-libs`传入兼容GNU Fortran的LP64/OpenMP链接选项。

编译日志、`link.log`和`build-report.json`保存在工作目录。修改构建设置时，使用独立工作目录。

## 3. 配置运行环境

设置许可证位置及外部运行库搜索路径：

```bash
export PANGULU_LICENSE_PATHS=/path/to/your/authorized-license.key
export LD_LIBRARY_PATH=/path/to/external/runtime-libraries:$LD_LIBRARY_PATH
```

可执行文件包含编译时选定的PanguLU和CUDA本机库搜索路径。BLAS、Fortran等运行库可能需要额外配置路径。

## 4. 运行小网格案例

```bash
python3 -m unittest discover -s tests -v
mkdir -p _work/runs
python3 tools/run_case.py dam --profile small \
  --binary _work/build/ccx_pangulu \
  --output _work/runs/dam-small
```

`small`模式使用上游网格，启用温度、热流及平衡输出。替换`dam`即可检查其他案例的小网格模型。

## 5. 运行提供的加密模型

```bash
python3 tools/run_case.py dam \
  --binary _work/build/ccx_pangulu \
  --output _work/runs/dam-full
```

可选案例名为`dam`、`ax6ht`、`capacitor`、`oneel20fi`和`pan`。运行工具解压输入模型后启动CalculiX。每次运行使用一个新输出目录。

请独占GPU并串行运行案例。组装与结果处理使用4线程，适配器后端使用40个主机线程。`small`模式的应用限时为240秒，加密模型为600秒，终止宽限均为10秒。

### 输出文件

| 文件 | 内容 |
|---|---|
| `model.inp` | CalculiX实际使用的输入模型 |
| `model.frd`、`model.dat`、`model.sta` | 原生物理场、表格结果及增量信息 |
| `application.stdout`、`application.stderr` | 应用消息与求解计时 |
| `application-exit.json` | 退出状态与应用运行时间 |
| `resources.txt` | GNU time的CPU、内存统计，可用时生成 |
| `actual.call*.csr.gz` | 每次调用的完整矩阵与右端项 |
| `log-validation.json` | 规模、调用次数、残差及计时检查 |
| `field-validation.json` | 温度、热流、边界及平衡检查 |
| `validation.json` | 运行检查的汇总状态 |

如需与同网格参考场比较，可添加`--reference-frd /path/to/reference.frd.gz`。计时定义与性能统计见[性能与验证](RESULTS.zh-CN.md)。

## 6. 可选：生成其他网格

提供的加密模型可直接运行。需要从上游案例生成网格时，使用：

```bash
python3 tools/generate_mesh.py dam --factor 32 --output _work/mesh-dam
python3 tools/generate_mesh.py ax6ht --factor 125 --output _work/mesh-ax6ht
python3 tools/generate_mesh.py capacitor --factors 3 16 400 --output _work/mesh-capacitor
python3 tools/generate_mesh.py oneel20fi --factor 44 --output _work/mesh-oneel20fi
python3 tools/generate_mesh.py pan --factor 16 --output _work/mesh-pan
```

工具在指定目录生成`model.inp`与网格元数据。细分保留父单元几何，并将集合和边界面映射到对应子单元。修改后的网格属于独立测试配置。

## 7. 可选：检查系数更新

`tests/changed_values.cpp`求解两个图结构相同、系数和右端项不同的3×3系统，两者的已知解均为(1,2,3)。

```bash
c++ -O2 -std=c++17 -fopenmp tests/changed_values.cpp \
  _work/build/obj/ccx_pangulu.cpp.o \
  -L/path/to/external/pangulu/lib -lpangulu -lpangulu_platform \
  -L/path/to/external/cuda/lib64 -l:libcudart.so.12 -lz -ldl \
  -o _work/changed_values
mkdir _work/reuse-check
HULL_CANONICAL="$PWD/_work/reuse-check/canonical.csr.gz" \
HULL_ACTUAL_PREFIX="$PWD/_work/reuse-check/actual" \
HULL_CREATE_CANONICAL=1 _work/changed_values
```

解满足检查条件时，输出会显示一次分析、两次数值分解及`CHANGED_VALUES_PASSED`。
