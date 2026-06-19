# BRU performance maps

Reconstructed **compressor** and **turbine** performance maps for the NASA **Brayton Rotating
Unit (BRU)** — a single-shaft turbomachine (4.25-in centrifugal compressor + 4.97-in radial-inflow
turbine + alternator) for a recuperated closed Brayton-cycle space-power system. The maps are built
from digitized 1969–1972 NASA rig reports and are formatted for an **NPSS-style off-design cycle
model** (they drop directly into the sibling [`../cycle`](../cycle) GasCycle solver).

Both maps are complete, validated, and integrated. This document is the entry point: what the
deliverables are, how to use them, and how to regenerate them.

---

## Start here — the two files you want

For a cycle model running at **nominal design conditions** (design fluid, pressure, temperature, and
Reynolds number), use exactly these two — and nothing else:

| Component | File | What it is |
|---|---|---|
| **Compressor** | **`maps/compressor/compressor_argon.map`** | clean aerodynamic R-line map |
| **Turbine** | **`maps/turbine/turbine_argon_tt.map`** | clean aerodynamic map, total-to-total PR (pairs with a cycle turbine element) |

Each is **self-contained**: the Reynolds correction is built in as a nested `S_Re` subelement, a no-op
at the **argon-rig reference** Reynolds number (`RNI = 1`) and active off it. The real He-Xe ~10 kWe
operational point sits at `RNI ≈ 1.6` (higher Re than the rig — the rig was set at the low end of the
He-Xe operating range), where the socket adds a small **favorable** ~+0.5 η-pt; near design it's
negligible. So these two files give you the clean rig aerodynamics — query at the design point and you
get compressor PR 1.91 / η 0.80 and turbine PR_tt 1.65 / η 0.92 (§"Validation") — and the socket nudges
efficiency to your operating Reynolds when you set `RNI`. Scale them to your cycle's design point and
pick the gas scalar as in "How to use the maps" below; that's the whole workflow.

Everything else under `maps/` is **situational** and *not* needed for design-point modeling — see the
table below. (Notably the `calibrated/indicated_*` maps deliberately bake in a hot-loop measurement
artifact and must **not** be used as aerodynamic maps.)

## All deliverables

| Need | File |
|---|---|
| **Compressor**, cycle modeling | `maps/compressor/compressor_argon.map` ← **the one** |
| **Turbine**, total-to-total cycle element | `maps/turbine/turbine_argon_tt.map` ← **the one** |
| Turbine, native total-to-static PR axis | `maps/turbine/turbine_argon.map` |
| Reynolds correction | built into the base maps (embedded `S_Re`); table dump in `*_reynolds_correction.csv` |
| Real-loop compressor efficiency band (situational) | `maps/compressor/calibrated/` (aero / indicated) |
| Reynolds / clearance scenarios, pre-baked (situational) | `maps/compressor/variants/`, `maps/turbine/variants/` |
| Portable long-format grids (non-NPSS tools) | `*_map_gridded.csv` next to each `.map` |

All maps are in **argon-equivalent corrected, normalized coordinates** and use the native NPSS
`.map` layout — a `Subelement <Type> S_map { … }` block carrying its design-scalar declarations and
its tables. The **compressor** is a `CompressorRlineMap` R-line map: tables `TB_Wc`, `TB_PR`, `TB_eff`
(leaf arrays `WcMap`, `PRmap`, `effMap`) as functions of `(alphaMap, NcMap, RlineMap)` with `RlineMap`
1.0 (surge) … 3.0 (max flow), design near `RlineMap ≈ 2.21`, two identical `alphaMap = 0/1` sub-maps
(fixed geometry — any alpha returns the same map), and the surge line declared as `RlineStall = 1.0`
(so NPSS computes the SMW/SMN stall margins). Each table carries its own per-axis interp/extrap methods
(`lagrange3` in speed/R-line/PR, `linear` in alpha and all extrapolation). The Reynolds correction is a
nested `Subelement {Compressor,Turbine}ReynoldsEffects S_Re { … }` embedded **inside** each base map's
`S_map` block (filling its `S_Re` socket; a `*_reynolds_correction.csv` dump sits alongside for
inspection). The **turbine** is a `TurbinePRmap` map: tables `TB_Wp`, `TB_eff` (leaf arrays `WpMap`,
`effMap`) as functions of `(NpMap, PRmap)` — pressure ratio is an independent, not an output.
`NcMap = N/N_des` (turbine `NpMap`, 100 % = 1.0), `WcMap = We/We_des` (turbine `WpMap`). See
`maps/README.md` for the folder layout.

## How to use the maps

1. **Load.** The `.map` files use NPSS `Table` syntax. GasCycle parses them with `read_npss_map`;
   the compressor R-line map converts to a `(Nc, Wc) → (PR, eff)` performance map, while the turbine
   map is read directly as `(Nc, PR) → (Wc, eff)`. The integration tests in `scripts/test/` are
   minimal worked examples.
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
     (both embed the Reynolds S_Re subelement via generate_reynolds_maps.{compressor,turbine}_s_re_block)
2. generate_map_variants.py      + reynolds_correction.py   -> maps/compressor/variants/
   generate_calibrated_variants.py + calibration.py         -> maps/compressor/calibrated/
   generate_turbine_variants.py  + turbine_corrections.py   -> maps/turbine/variants/
   generate_reynolds_maps.py     + reynolds/turbine_corr.   -> maps/*/{,turbine_}reynolds_correction.csv
3. test/test_{map,turbine}_integration.jl   load the maps into GasCycle (needs ../cycle + Julia)
```
`corrected_params.py` holds shared design constants; `reynolds_correction.py`, `calibration.py`,
`turbine_corrections.py` are imported libraries, not run directly. `generate_reynolds_maps.py` is both
importable (the stage-1 builders call it to embed the `S_Re` socket) and runnable (writes the
`*_reynolds_correction.csv` dumps). To rebuild everything: run stage 0, then the two stage-1 builders,
then the three stage-2 generators.

**Source data.** The digitized inputs in `data/digitized/` come from the user's manual marker-center
WebPlotDigitizer reads (`data/manual_digitize/`), which are the validated ground truth — cross-checked
three independent ways and against each figure's printed design dot. `data/raw_text/` holds extracted
report text. The 23 source PDFs live one level up in the parent corpus.

## Validation summary & known limits

- **Compressor** queries at design to PR 1.912 / η 0.798 (spec 1.90 / 0.795 — ~1 % residual from the
  Table-II design flow sitting just right of the plotted marker). Validated vs hot closed-loop
  TM X-2350 (CRP 80 %, in-loop 74.5 %, rotor power balance closes to 0.67 %). Surge line
  (`RlineStall = 1.0`) matches the digitized Fig-8 surge curve to RMS ≈ 0.036 / max ≈ 0.062 in PR.
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
