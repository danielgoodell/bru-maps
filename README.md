# BRU performance maps

Reconstructed **compressor** and **turbine** performance maps for the NASA **Brayton Rotating
Unit (BRU)** — a single-shaft turbomachine (4.25-in centrifugal compressor + 4.97-in radial-inflow
turbine + alternator) for a recuperated closed Brayton-cycle space-power system. The maps are built
from digitized 1969–1972 NASA rig reports and are formatted for an **NPSS-style off-design cycle
model** (they drop directly into the sibling [`../cycle`](../cycle) GasCycle solver).

Both maps are complete, validated, and integrated. This document is the entry point: what the
deliverables are, how to use them, and how to regenerate them.

---

## Deliverables — what to use

| Need | File |
|---|---|
| **Compressor**, cycle modeling | `maps/compressor/compressor_argon.map` |
| **Turbine**, total-to-total cycle element (recommended) | `maps/turbine/turbine_argon_tt.map` |
| Turbine, native total-to-static PR axis | `maps/turbine/turbine_argon.map` |
| Real-loop compressor efficiency band | `maps/compressor/calibrated/` (aero / indicated) |
| Reynolds / clearance scenarios | `maps/compressor/variants/`, `maps/turbine/variants/` |
| Portable long-format grids (non-NPSS tools) | `*_map_gridded.csv` next to each `.map` |

All maps are in **argon-equivalent corrected, normalized coordinates**:
`NcMap = N/N_des`, `WcMap = We/We_des`, `betaMap` = R-line (0 = surge/low-PR … 1 = max-flow/high-PR),
giving `PRmap`, `effMap`, `WcMap` tables. See `maps/README.md` for the folder layout.

## How to use the maps

1. **Load.** The `.map` files use NPSS `Table` syntax. GasCycle parses them via
   `read_npss_map → to_performance_map → query(Nc, Wc)`; the integration tests in `scripts/test/`
   are minimal worked examples.
2. **Scale to your design point.** The maps are normalized, so set the element's map scalers
   (`s_Nc, s_Wc, s_PR, s_eff`) so that at your cycle design point the corrected speed/flow/PR/eff
   equal the BRU design values:

   | | Compressor | Turbine |
   |---|---|---|
   | Equivalent flow We_des | 0.581 lb/s | 0.486 lb/s |
   | Equivalent speed N_des | 51 176 rpm | (design ν = 0.70) |
   | Design PR | 1.90 (P6/P1) | 1.69 total-to-static / 1.645 total-to-total |
   | Design efficiency | 0.795 (peak 0.819) | 0.913 (total, 1→3) |

   These constants live in `scripts/build/corrected_params.py`.
3. **Pick the right gas.** Because all candidate fluids are monatomic (γ = 5/3), the **normalized
   map is gas-independent** — it serves argon, He-Xe (MW 83.8), and krypton. The fluid enters only
   through the per-gas design scalars (flow/speed scale as √MW). One turbine caveat: its ν-from-PR
   formula bakes in argon PR_des = 1.69; for He-Xe compute ν from actual conditions. See the
   gas-transfer sections of `docs/compressor.md` / `docs/turbine.md`.
4. **Choose a fidelity level.** Use the base `.map` for clean aerodynamic off-design. For real-loop
   representation, pick a `calibrated/` efficiency band (compressor) and/or a `variants/` member
   (Reynolds, clearance). **Heat leak is modeled at the cycle level** (Option B: clean aero maps +
   a turbine→compressor heat term in the cycle), not baked into the maps.

## How it's built (regenerate)

The pipeline is in `scripts/build/` (run from there), with solver integration tests in
`scripts/test/`. Stages:

```
0. ingest_manual_digitize.py     data/manual_digitize/ (user WPD reads) -> data/digitized/*.csv
1. build_compressor_map.py       data/digitized/comp_*      -> maps/compressor/compressor_argon.map
   build_turbine_map.py          data/digitized/turbine_*   -> maps/turbine/turbine_argon{,_tt}.map
2. generate_map_variants.py      + reynolds_correction.py   -> maps/compressor/variants/
   generate_calibrated_variants.py + calibration.py         -> maps/compressor/calibrated/
   generate_turbine_variants.py  + turbine_corrections.py   -> maps/turbine/variants/
3. test/test_{map,turbine}_integration.jl   load the maps into GasCycle (needs ../cycle + Julia)
```
`corrected_params.py` holds shared design constants; `reynolds_correction.py`, `calibration.py`,
`turbine_corrections.py` are imported libraries, not run directly. To rebuild everything: run stage 0,
then the two stage-1 builders, then the three stage-2 generators.

**Source data.** The digitized inputs in `data/digitized/` come from the user's manual marker-center
WebPlotDigitizer reads (`data/manual_digitize/`), which are the validated ground truth — cross-checked
three independent ways and against each figure's printed design dot. `data/raw_text/` holds extracted
report text. The 23 source PDFs live one level up in the parent corpus.

## Validation summary & known limits

- **Compressor** queries at design to PR 1.912 / η 0.798 (spec 1.90 / 0.795 — ~1 % residual from the
  Table-II design flow sitting just right of the plotted marker). Validated vs hot closed-loop
  TM X-2350 (CRP 80 %, in-loop 74.5 %, rotor power balance closes to 0.67 %).
- **Turbine** queries at design to PR_ts 1.694 / η 0.916 (spec 1.69 / 0.913); design flow now
  coincides with the 100 % line (−0.2 % vs Table I, the old figure-vs-table offset resolved by the
  manual reads). Cross-validated vs in-loop TM X-2350 (clearance) and He-Xe TM X-67996 (gas transfer).
- **Gas transfer** to He-Xe / Kr validated for both components (√MW scaling; normalized map invariant).
- **Limits:** efficiency carries ~±0.01–0.02; low-Re surge/max-flow shift is not applied (reserved
  for loop calibration); the cycle-level turbine→compressor heat-leak term is a modeling choice left
  to the cycle deck.

## Documentation

- `docs/compressor.md`, `docs/turbine.md` — per-component: design point, corrected parameters, map
  construction, corrections, gas transfer, validation, integration.
- `docs/inventory.md` — the 23 source reports, triaged by role.
- `maps/README.md` — maps folder layout and which file to use.
