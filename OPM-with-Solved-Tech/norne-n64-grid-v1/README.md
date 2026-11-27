# Norne N64

Norne N64 is a 4 × 4 × 4 refinement of the public Norne model. Its 184 × 448 × 88 grid has 2,843,196 active cells. The Eclipse input preserves the historical schedule from 1997-11-06 to 2006-12-01 (3,312 simulated days).

## Downloads

| File | Contents | Download size |
| --- | --- | ---: |
| [Norne_N64_4x4x4.DATA.tar.gz](https://raw.githubusercontent.com/solved-tech-dev/case-studies/main/OPM-with-Solved-Tech/norne-n64-grid-v1/Norne_N64_4x4x4.DATA.tar.gz) | NORNE_N64.DATA and GRID_PROPS.INC | 99.9 MiB |
| [Norne_N64_4x4x4.EGRID.gz](https://raw.githubusercontent.com/solved-tech-dev/case-studies/main/OPM-with-Solved-Tech/norne-n64-grid-v1/Norne_N64_4x4x4.EGRID.gz) | Eclipse binary corner-point grid | 44.8 MiB |

Keep the two simulator input files in the same directory after extraction.

```sh
tar -xzf Norne_N64_4x4x4.DATA.tar.gz
flow Norne_N64/NORNE_N64.DATA
```

To extract the binary grid, run `gzip -dk Norne_N64_4x4x4.EGRID.gz`. The extracted grid is 265,494,912 bytes.

## SHA-256

| File | SHA-256 |
| --- | --- |
| `Norne_N64_4x4x4.DATA.tar.gz` | `777ca3615d9c7151ba81dd0e95c4e72abe6b9c33463c39815723eb312354772f` |
| `Norne_N64_4x4x4.EGRID.gz` | `9b726288a9ecf5626a176036bc6f1775f028721c0d1eb221a4764869acdec735` |
| `Norne_N64_4x4x4.EGRID` | `f11a4dff73814028198dbb65a2728854c9b465eecb808bc598ce5135ccef80bb` |

See the [Norne N64 case notes](../examples/norne/README.md) for model details and shorter runs.

Derived from the [public Norne ATW 2013 deck](https://github.com/OPM/opm-data/tree/master/norne), copyright Statoil (2015). See the upstream [Open Database License](https://opendatacommons.org/licenses/odbl/1-0/) and [Database Contents License](https://opendatacommons.org/licenses/dbcl/1-0/).
