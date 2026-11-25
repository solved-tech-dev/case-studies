# Sources and licensing · 来源与许可

The Solved Tech adapter, scripts and documentation follow the licensing of the other with-Solved-Tech application packages, **AGPL-3.0-only**. See [LICENSE](LICENSE). The upstream portions of the patch and the adapted FrontISTR tutorial inputs retain the [FrontISTR Commons MIT notice](licenses/FrontISTR-MIT.txt).

FrontISTR 5.9 source revision: [8ad8993860151e7070f0ca6078b8a0157966561e](https://github.com/FrontISTR/FrontISTR/tree/8ad8993860151e7070f0ca6078b8a0157966561e). Model sources: [19_conrod](https://github.com/FrontISTR/FrontISTR/tree/8ad8993860151e7070f0ca6078b8a0157966561e/tutorial/19_conrod), [01_elastic_hinge](https://github.com/FrontISTR/FrontISTR/tree/8ad8993860151e7070f0ca6078b8a0157966561e/tutorial/01_elastic_hinge), [16_heat_block](https://github.com/FrontISTR/FrontISTR/tree/8ad8993860151e7070f0ca6078b8a0157966561e/tutorial/16_heat_block). Each input directory contains a refined mesh and a direct-solver control file. The two tetrahedral meshes were refined by edge bisection. Each original hexahedron of the heat model was split into 27 cells. The figures visualize nodal solutions of the published meshes.

解得科技新增的适配代码、脚本和文档沿用同系列项目的 AGPL-3.0-only 许可。三个模型来自上述固定版本的官方教程。网格经加密，控制文件选择直接求解器。上游代码和教程保留原 MIT 版权声明。

PanguLU 26.48, its license, the CUDA toolkit and third-party dependencies are obtained separately and are subject to their respective terms. [PanguLU SDK and licensing / 软件与授权](https://solved-tech.com/)
