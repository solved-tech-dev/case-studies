# 构建与运行

## 环境要求

准备与验证脚本需要 Python 3.10 或以上版本，依赖见 `requirements.txt`。

使用 Linux x86_64、NVIDIA GPU、兼容 CUDA 12.9 的驱动与完整 CUDA 工具包，以及 **PanguLU 26.46 ABI 2**。接口按 2026-11-12 更新的 26.46 包整理，原包 SHA-256 为 `6aeb44706e01827e56924e2ea4d97068276ddbf50418fdf03103f69582dfff29`。SDK 需单独取得，本文不嵌入带授权信息的下载链接。

MOOSE 需要实数双精度 PETSc 及相匹配的 MPI、libMesh。开发核验采用 PETSc 3.25.4；插件必须与应用实际链接的 PETSc 一致。编译器和动态加载器需能找到 CUDA 头文件、运行库、cuBLAS 及 SDK 所需依赖。设备许可证另行获取，放置到 `~/.pangulu/license.key`，仅允许本人访问；不要放入本仓库。

MOOSE 依赖安装参见 [官方安装说明](https://mooseframework.inl.gov/getting_started/installation/index.html)。脚本不会改动现有安装或自动升级 PETSc、MPI、CUDA。

## 1. 获取上游并构建应用

以下命令在本代码包根目录执行。`/path/to/...` 替换为自己的安装位置。

```bash
python3 -m pip install -r requirements.txt
python3 tools/fetch_upstream.py _upstream/moose
export MOOSE_DIR="$PWD/_upstream/moose"
# 按 MOOSE 源码安装说明构建兼容的 PETSc 和 libMesh，随后指定其位置：
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-your-build  # 前缀安装可留空
export LIBMESH_DIR=/path/to/libmesh/installed
make -C app -j 8 MOOSE_DIR="$MOOSE_DIR" METHOD=opt
```

上游固定为 `07b49de938ed357c136b487078065406d838aa26`。获取脚本只创建新的检出目录，不覆盖已有源码。依赖及必要子模块需按正常 MOOSE 源码构建流程准备。`app/` 只注册原有 Solid Mechanics、Heat Transfer 等所需模块，不替换物理实现。已有的实数双精度 MOOSE 应用若注册了这些模块且支持相同输入，也可直接使用。

## 2. 构建 PanguLU 插件

```bash
export PANGULU_ROOT=/path/to/pangulu-26.46
export CUDA_HOME=/path/to/cuda-12.9
export PATH="$CUDA_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$PANGULU_ROOT/lib:$CUDA_HOME/lib64:${LD_LIBRARY_PATH:-}"
python3 tools/build_plugin.py
```

默认编译器为 `mpicxx`，也可通过 `--cxx` 指定与 PETSc 相同的 MPI 编译器。`--petsc-dir`、`--petsc-arch`、`--pangulu-root`、`--cuda-root` 可替代环境变量。生成文件为 `build/libmoose_pangulu.so`。若 CUDA 分散在多个软件包目录中，需通过编译器封装和动态库搜索路径补齐头文件与库。

不要混用不同安装的 PETSc 头文件和库；运行前可用 `ldd build/libmoose_pangulu.so` 确认没有缺失依赖。

## 3. 运行完整应用

在已分配的 GPU 资源上执行。当前为单 MPI 进程，命令不需要 `mpirun`。

```bash
export OMP_NUM_THREADS=8
export OMP_DYNAMIC=FALSE
export OPENBLAS_NUM_THREADS=1
python3 tools/run_case.py thermomech \
  --exe "$PWD/app/moosesolved-opt" \
  --profile small --output runs/thermomech-small
```

输出目录必须是新目录，以免覆盖旧结果。程序准备必要的网格，运行完整物理模型，并检查收敛、物理时间、完整字段及插件输出的独立原 CSR 残差。保存 `run.json`、网格和应用命令、`application.log`、Exodus 物理场、CSV 与 `validation.json`。

改为 `--profile million` 即使用已提供的百万阶网格。大网格比小规模核验需要更多内存，请据环境分配资源。运行器不绑定集群、硬件型号或固定调度队列，也不包含预设性能目标或自动运行超时。

六个案例的统一入口：

```bash
for case in thermomech annular-torsion mug pinched-cylinder gap-heat hex20-heat; do
  python3 tools/run_case.py "$case" --exe "$PWD/app/moosesolved-opt" \
    --profile million --output "runs/$case-million"
done
```

这会按顺序各运行一次完整应用。`--prepare-only` 仅生成网格和 `run.json`，不求解物理问题；网格格式转换仍会调用 MOOSE。

## 4. 参考结果与测试

小网格可以使用 PETSc 自带顺序 LU 作为参考：

```bash
python3 tools/run_case.py thermomech --exe "$PWD/app/moosesolved-opt" \
  --profile small --solver reference --reference-solver petsc \
  --output runs/thermomech-reference
python3 tools/validate.py runs/thermomech-small --reference runs/thermomech-reference
c++ -std=c++17 -O2 -Iadapter tests/test_sparse.cpp -o /tmp/moose-test-sparse
/tmp/moose-test-sparse
python3 tools/build_plugin.py --native-test --output build/test-native
./build/test-native
python3 tools/build_plugin.py --petsc-test --output build/test-petsc
./build/test-petsc -dll_prepend "$PWD/build/libmoose_pangulu.so"
```

原生测试需有许可证的 GPU，检查显式零与非零之间的值更新、同矩阵分解复用、图变化后的重分析及真正的两列 RHS 块求解。测试不比较运行时间。内置 `petsc` LU 参考适用于小网格；也可明确选择安装中已有的其他 PETSc LU 解法器。

详见[结果验证](VALIDATION.zh-CN.md)及[接入说明](INTEGRATION.zh-CN.md)。
