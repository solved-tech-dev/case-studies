# FEBio-with-Solved-Tech

[简体中文](README.zh-CN.md) · [Build and run](docs/BUILD-AND-RUN.md) · [Solver interface](docs/INTEGRATION.md)

Run three refined FEBio verification problems with PanguLU 26.46. FEBio assembles the tangent matrix, advances the load and time steps, and writes the native plot file. A loadable FECore plugin handles each linear system. The FEBio source tree needs no edits.

| Case | Model | HEX8 elements | Equations |
| --- | --- | ---: | ---: |
| [bp10](cases/bp10/README.md) | Pressurized biphasic quarter-cylinder | 276,480 | 904,176 |
| [tu01](cases/tu01/README.md) | Uniaxial extension of a hyperelastic block | 256,000 | 803,518 |
| [mi06](cases/mi06/README.md) | In-plane loading of a neo-Hookean plate | 202,500 | 1,017,000 |

Start with the build guide to generate a model from the pinned FEBio TestSuite, compile the plugin and open the resulting displacement and stress fields in FEBio Studio. The preparation script verifies each generated input against the exact file used in the application runs. [Dependencies and licensing](NOTICE.md)
