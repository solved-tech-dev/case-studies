# Norne N64 Benchmark Case Notes

Norne N64 is a 4 × 4 × 4 global refinement of the public Norne field model. The logical grid grows from 46 × 112 × 22 to **184 × 448 × 88**, with **2,843,196 active cells** and **8,529,588 scalar unknowns** in the three-component black-oil reservoir system.

## Simulation input

Download [`Norne_N64_4x4x4.DATA.tar.gz`](https://raw.githubusercontent.com/solved-tech-dev/case-studies/main/OPM-with-Solved-Tech/norne-n64-grid-v1/Norne_N64_4x4x4.DATA.tar.gz) (99.9 MiB) from the [Norne N64 data directory](../../norne-n64-grid-v1/README.md). The deck contains the historical schedule from **1997-11-06** to **2006-12-01** (3,312 simulated days).

```sh
tar -xzf Norne_N64_4x4x4.DATA.tar.gz
flow Norne_N64/NORNE_N64.DATA
```

The archive contains these two files in the same directory:

| File | Contents | Uncompressed size |
| --- | --- | ---: |
| `NORNE_N64.DATA` | Eclipse input with fluid tables, well connections, controls and dates | 2.48 MiB |
| `GRID_PROPS.INC` | N64 corner-point grid, active-cell mask and rock properties | 720.22 MiB |

Keep both files together so that `INCLUDE` resolves from the main deck. For a shorter benchmark, use a copy of the main deck and end its `SCHEDULE` at the desired reporting date. The same `GRID_PROPS.INC` can be reused.

To use the solver integrations in this repository, pass the path to `NORNE_N64.DATA` as `--deck` to `scripts/run-case.sh` with a configuration from `examples/configs/`.

Archive SHA-256: `777ca3615d9c7151ba81dd0e95c4e72abe6b9c33463c39815723eb312354772f`.

## Corner-point Grid Files

- **Original Norne grid (46 × 112 × 22).** The upstream [`NORNE_ATW2013.DATA`](https://github.com/OPM/opm-data/blob/master/norne/NORNE_ATW2013.DATA) reads the geometry from [`IRAP_1005.GRDECL`](https://github.com/OPM/opm-data/blob/master/norne/INCLUDE/GRID/IRAP_1005.GRDECL) and the active-cell mask from [`ACTNUM_0704.prop`](https://github.com/OPM/opm-data/blob/master/norne/INCLUDE/GRID/ACTNUM_0704.prop).
- **Refined Norne N64 grid (184 × 448 × 88).** Download [`Norne_N64_4x4x4.EGRID.gz`](https://raw.githubusercontent.com/solved-tech-dev/case-studies/main/OPM-with-Solved-Tech/norne-n64-grid-v1/Norne_N64_4x4x4.EGRID.gz) (44.8 MiB). Extract with `gzip -dk Norne_N64_4x4x4.EGRID.gz` to obtain the Eclipse binary grid containing `COORD`, `ZCORN` and `ACTNUM`. Uncompressed SHA-256: `f11a4dff73814028198dbb65a2728854c9b465eecb808bc598ce5135ccef80bb`.

The N64 model is derived from the public Norne ATW 2013 deck, copyright Statoil (2015). The upstream notices specify the [Open Database License](https://opendatacommons.org/licenses/odbl/1-0/) and the [Database Contents License](https://opendatacommons.org/licenses/dbcl/1-0/).
