# FrontISTR 接入 PanguLU

[English](README.md) · [构建与运行](docs/BUILD.md) · [来源与许可](NOTICE.md)

用 FrontISTR 5.9 组装真实有限元方程，在单张 CUDA GPU 上调用 PanguLU 26.48 求解，再由 FrontISTR 输出位移、应力或温度场。包内附带推送中三个案例的完整网格、控制文件和结果校验值。

| 案例 | 节点 / 单元 | 方程数 | 计算内容 |
|---|---:|---:|---|
| [三维连杆](cases/conrod) | 333,340 / 226,614 | 1,000,020 | 大端固定，小端施加位移，计算位移和节点应力。[位移图](cases/conrod/field.png) |
| [弹性铰链](cases/hinge) | 333,340 / 227,819 | 1,000,020 | 固定两组边界，沿第一坐标方向施加节点载荷。[位移图](cases/hinge/field.png) |
| [块体热传导](cases/heat) | 914,550 / 868,320 | 914,550 | 中心与四组侧面分别给定 100 和 20 的温度，求解温度相关材料的热传导问题。[温度剖面](cases/heat/field.png) |

**接入位置。** [上游补丁](patches/frontistr-pangulu.patch) 在 HEC-MW 的 MUMPS 入口增加可选分支。控制文件仍使用 `METHOD=MUMPS`；运行脚本设置 `FRONTISTR_SOLVER=pangulu` 后，求解阶段转入 [bridge.cpp](adapter/bridge.cpp)。按单进程、单 GPU、实数双精度、单右端项工作。未启用该环境变量时，走原版 MUMPS 分支。

**复用与校验。** 适配器展开对称存储、合并重复项，再将 CSR 转为 CSC 送入 [PanguLU 后端](adapter/gpu_backend.inc)。结构相同就复用分析，系数也相同就复用数值分解。热传导案例的六次求解做一次分析、一次初始分解和五次重分解。每次按原方程独立计算相对残差。运行脚本还核对调用次数、复用记录及最终场文件的节点数、各字段和参考范数。

推送使用 PanguLU 26.48（CUDA 12.9）和一张 A100 PCIe 40 GB。示例默认以 4 个 OpenMP 线程运行 FrontISTR，每次进入 PanguLU 求解阶段时切换至 40 个主机线程，返回时恢复。可通过 `FRONTISTR_PANGULU_THREADS` 调整求解线程数。

构建后可直接执行下列命令。每次使用新的输出目录，日志和物理场都保存在该目录内。

```bash
python3 scripts/run.py conrod runs/conrod --fistr1 build/frontistr/fistr1/fistr1
python3 scripts/run.py hinge  runs/hinge  --fistr1 build/frontistr/fistr1/fistr1
python3 scripts/run.py heat   runs/heat   --fistr1 build/frontistr/fistr1/fistr1
```
