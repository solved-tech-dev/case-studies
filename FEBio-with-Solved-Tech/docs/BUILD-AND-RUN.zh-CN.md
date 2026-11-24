# 生成网格、构建插件与运行案例

[English](BUILD-AND-RUN.md) · [接口和复用规则](INTEGRATION.zh-CN.md)

以下命令从 `FEBio-with-Solved-Tech/` 根目录执行。生成网格需要 Python 3 和 NumPy。构建、运行插件需要 Linux x86-64、支持 C++17 与 OpenMP 的 GCC、匹配的 FEBio 4.13.0 头文件与共享库、PanguLU 26.46（ABI 2，CUDA 12.9 版）、对应 CUDA 工具链，以及可供 PanguLU 使用的 GPU 和授权。集成验证使用单张 A100 40 GB；实测插件在求解器调用期间将 OpenMP 线程数设为 40，返回 FEBio 时恢复原设置。

## 1. 生成输入

从上游固定版本获取三个原始模型。生成程序按原形函数或网格坐标细分父单元、映射边界条件，并将求解器入口设为 `febio_gpu_lu`。程序先核对原始输入，再核对输出文件与实测输入的 SHA-256。可在命令末尾只指定一个案例。

```bash
git clone --filter=blob:none --no-checkout https://github.com/febiosoftware/TestSuite.git TestSuite
git -C TestSuite sparse-checkout set --no-cone /Verify/bp10.feb /Verify/tu01.feb /Verify/mi06.feb
git -C TestSuite checkout --detach 986005fc1865f5de7f2c789b8153d6f8fe2c0acd
python3 -m pip install numpy
python3 scripts/prepare.py TestSuite bp10 tu01 mi06
```

## 2. 编译插件

下面是示例安装路径，请按本机安装位置替换。FEBio 可执行文件须与用于编译插件的 FECore 兼容。先按供应商说明完成 PanguLU 的运行环境与授权配置。

```bash
export FEBIO_SRC=/opt/FEBio/include
export FEBIO_LIB=/opt/FEBio/lib
export FEBIO_BIN=/opt/FEBio/bin/febio4
export PANGULU_ROOT=/opt/pangulu-26.46
export CUDA_HOME=/usr/local/cuda-12.9
export LD_LIBRARY_PATH="$FEBIO_LIB:$PANGULU_ROOT/lib:$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
bash scripts/build.sh
```

`FEBIO_SRC` 指向包含 `FECore/` 的目录。若 CUDA 的运行时头文件、编译器头文件和库分开安装，可分别设置 `CUDA_INCLUDE`、`CUDA_NVCC_INCLUDE` 和 `CUDA_LIB`。编译产物为 `FEBio-with-Solved-Tech/build/pangulu.so`。

## 3. 运行和检查

每个案例使用一个新的输出目录。保留 `console.log`、`application.log`、`physical.xplt`。脚本检查 FEBio 正常终止、存在原方程残差记录且生成了非空的物理场文件。

```bash
mkdir -p runs
bash scripts/run.sh bp10 runs/bp10
bash scripts/run.sh tu01 runs/tu01
bash scripts/run.sh mi06 runs/mi06
```

`console.log` 中的 `FEBIO_RESIDUAL` 对应每一次实际求解。插件检查 ‖Ax−b‖₂/‖b‖₂ ≤ 10⁻⁸；当 b 为零时使用绝对残差。用 FEBio Studio 打开 `physical.xplt` 可查看各收敛时刻的位移、应力等场。

如已用另一求解器生成同一网格、材料、载荷及输出设置下的参考 XPLT，可逐场核对。例如：

```bash
python3 scripts/compare_fields.py cases/tu01/tu01-8.feb \
    /path/to/reference/physical.xplt runs/tu01/physical.xplt
```

对照脚本读取本组案例使用的 FEBio 未压缩 XPLT 3.5，检查节点、单元、输出变量和收敛时刻；所有场分量均要求 |候选值−参考值| / max(1, |参考值|) ≤ 10⁻⁶。不一致时返回非零退出码。脚本会将两份场数据读入内存。
