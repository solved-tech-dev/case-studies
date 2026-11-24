# 构建与运行

使用 Linux x86-64 和 Python 3.11 或更新版本。构建 ngspice 需要 Git、C 编译器、GNU Make、Autoconf、Automake、Libtool、Flex 和 Bison。Debian 或 Ubuntu 对应的构建包包括 `build-essential git autoconf automake libtool flex bison`。PanguLU SDK 和 CUDA 仅在回放阶段使用。

命令从代码包根目录执行。每次运行使用新的输出目录，既有结果可以继续保留。

## 运行 ngspice 并采集原方程

```bash
JOBS=8 ./tools/build-ngspice.sh build
case_id=pdn-interleaved-segments-dc
python3 tools/generate_cases.py "$case_id" --output-dir netlists
python3 tools/run_case.py capture "$case_id" \
    --ngspice build/ngspice-build/src/ngspice \
    --netlist "netlists/$case_id.cir" \
    --output-dir "runs/$case_id/klu"
```

构建脚本在独立目录检出 [upstream.lock.json](../upstream.lock.json) 固定的 ngspice 版本并加上采集接口。机器无法联网时，可把含该提交的本地 Git 仓库作为脚本的第二个参数。原生运行输出包括 `ngspice.log`、完整工作点结果 `field.raw`、`systems/` 中的原始方程和参考解，以及采集校验记录 `capture.json`。

已有启用 KLU 的标准 ngspice 时，可直接运行 `ngspice -n -b -r field.raw netlists/<case>.cir`。要导出方程，请使用本包构建的版本。

## 使用 PanguLU 回放

准备适用于 CUDA 12.9 的 PanguLU 26.46 ABI 2 SDK、运行库，以及所用 GPU 对应的有效许可证。按 SDK 的安装说明配置许可证。设置 SDK 和 CUDA 路径。CUDA 安装分设目录时，可以用 `CUDA_RUNTIME_INCLUDE`、`CUDA_RUNTIME_LIB_DIR` 和 `CUDA_COMPILER_INCLUDE` 指定对应位置。SDK 的库及其依赖目录按实际安装情况加入 `LD_LIBRARY_PATH`。

```bash
export PANGULU_ROOT=/path/to/pangulu-26.46
export CUDA_HOME=/path/to/cuda-12.9
./tools/build-pangulu.sh build
case_id=pdn-interleaved-segments-dc
CUDA_VISIBLE_DEVICES=0 python3 tools/run_case.py replay "$case_id" \
    --pangulu-replay build/pangulu-replay \
    --capture-dir "runs/$case_id/klu" \
    --output-dir "runs/$case_id/pangulu" --host-threads 1
```

原测试的 GPU 作业分配 40 个主机 CPU，PanguLU 的主机 OpenMP 线程数设为 1。程序使用当前可见的 CUDA 设备。`pangulu.log` 记录设备、同步后的求解器阶段和逐次校验。`result.json` 在三次求解全部校验通过后生成，分别报告求解器阶段累计耗时与回放进程总耗时。

在新的网表目录用 `python3 tools/generate_cases.py all --output-dir netlists` 生成五个完整网表。随后按[案例说明](../cases/README.zh-CN.md)中的 ID 逐一采集和回放。采集与回放可分开在 CPU 和 GPU 机器上进行，转移数据时请保留完整的采集目录。参考解核对记录见[原生方程采集验证](../validation/README.zh-CN.md)。
