# License and attribution

The unmodified files under `cases/*/upstream/` come from Idaho National Laboratory's MOOSE project. Their original paths, pinned revision and file digests are recorded in `upstream.lock.json`. They retain the upstream LGPL-2.1 license and copyright notices; see `licenses/MOOSE-LGPL-2.1.txt` and `licenses/MOOSE-COPYRIGHT`. The source repository is https://github.com/idaholab/moose.

Solved Tech's new adapter, application registration, mesh preparation scripts, case profiles, tests and documentation are provided under **AGPL-3.0-only**; see `LICENSE`. The application calls the original MOOSE physics objects. This package does not replace MOOSE's element kernels, materials or boundary-condition implementations.

MOOSE, PETSc, libMesh, MPI and their dependencies retain their own licenses. PanguLU and NVIDIA CUDA are separately licensed dependencies and are not redistributed here. Obtain a PanguLU SDK and a license for your GPU from [Solved Tech](https://solved-tech.com/). No SDK library, SDK header or device license is included in this source package.

# 许可与来源

`cases/*/upstream/` 中的原始输入与种子网格来自 MOOSE，保持原文件内容和上游许可。精确版本与文件来源见 `upstream.lock.json`，上游许可及版权说明见 `licenses/`。解得科技新增的适配、应用注册、网格处理、配置、测试与文档采用 AGPL-3.0-only。PanguLU SDK、CUDA 与设备许可证需另行获取，不随本代码包分发。
