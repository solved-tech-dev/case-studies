# 接入接口与输出

[English](INTEGRATION.md)

`tools/prepare.py` 在 PHG 求解器注册表加入 `pangulu`，并加入可选的原生场导出入口。PHG 和 PHG-Solid 的案例源码、材料、网格读取、组装以及边界处理来自官方归档。案例配置由 `tools/run_case.py` 写入工作目录。

`integration/phg/solver-pangulu.c` 通过 PHG 的求解器回调读取原矩阵的全部存储项，形成 CSR 和右端项。它调用 `adapter/pangulu_bridge.cpp` 中的 `phg_pangulu_solve_csr()`。适配器将 CSR 转为 CSC，保留全部存储项，用 `pangulu_execute()` 依次执行分析、数值分解和求解，再把解写回 PHG 的向量。三例使用实数双精度通用 LU，每次应用调用有一个右端项。自适应网格每轮重新分析和分解。

日志中的原生阶段记录用于诊断。运行报告检查各阶段返回码、矩阵尺寸、原 CSR 相对残差和完整物理场的有限性。相对残差定义为 `||Ax − b||₂ / ||b||₂`，阈值为 `1e-8`。固体案例输出 `Displacement`，球形界面输出 `u1_h` 和 `u2_h`。指定 `--capture-linear` 可另存每次调用的原 CSR、右端项及完整解。大网格的这些文件需要额外磁盘空间。

开发时可用 `tools/validate_connections.py --candidate pangulu --output _work/check-gpu` 在 GPU 上将三个小案例的完整场与本项目的 CPU 测试解对照。`--reference-only` 构建用于 CPU 接口检查。
