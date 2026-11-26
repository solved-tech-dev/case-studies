# Norne N64 Benchmark Case Notes

The Norne N64 case is a 4x4x4 global refinement of the public Norne field model:

- Base Norne grid: `46 x 112 x 22`
- Refined N64 logical grid: `184 x 448 x 88` (`7,254,016` logical cells)
- Active cells: `2,843,196`
- Black-oil 3-component main Jacobian (`BSR3`, FP64):
  - `2,843,196` block rows (`8,529,588` scalar unknowns)
  - `19,947,562` block nonzeros (`179,528,058` stored scalar entries in reservoir matrix)
  - Plus `36` standard well equations (`2,843,232` pressure system rows, `20,012,068` CSR nonzeros in CPRW)
- Simulation start date: `1997-11-06`
  - Day-1 window: `1997-11-06` to `1997-11-07` (`1.0` physical day)
  - 176-day window: `1997-11-06` to `1998-05-01` (`176.0` physical days)

## Obtaining the Data Deck

Large `.DATA` / `.INC` grid files and `.UNSMRY` / `.SMSPEC` / `.UNRST` outputs are intentionally excluded from this public repository.

To run Norne N64 locally:

1. Obtain the upstream Norne Eclipse deck from the [OPM opm-data repository](https://github.com/OPM/opm-data/tree/master/norne) (`NORNE_ATW2013.DATA`).
2. Apply 4x4x4 global refinement (or use your prepared refined `NORNE1.DATA` / `NORNE176.DATA` deck with shared `GRID_PROPS.INC`).
3. Pass the deck and one of the route configurations in `examples/configs/` to `scripts/run-case.sh`.

For day-1 validation without proprietary data, `reference-day1-summary.json` records the verified physical summaries, job IDs, and solver statistics from each route's reference runs.
