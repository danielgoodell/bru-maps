# Compressor map — 4.25-in sweptback centrifugal compressor

Primary source: **NASA TM X-2129 (1970)**, *Overall Performance in Argon of the 4.25-Inch
Sweptback-Bladed Centrifugal Compressor* (Tables I & II, Fig 8 PR-flow + surge, Fig 10 efficiency).

Deliverable: `maps/compressor/compressor_argon.map` (+ gridded CSV, validation PNG, and the
Reynolds / calibrated variant sets). Built by `scripts/build/build_compressor_map.py` from the
digitized inputs in `data/digitized/comp_*_argon.csv`.

---

## 1. Design point (two equivalent statements of one point)

| Quantity | He–Xe (design fluid, Table I) | Argon (rig fluid, Table II) |
|---|---|---|
| Working-fluid MW | 83.8 | 39.95 |
| Inlet total pressure P₁ | 13.5 psia | 20.25 psia |
| Inlet total temperature T₁ | 540 °R (300 K) | 540 °R |
| Weight flow W | 0.756 lb/s | 0.785 lb/s |
| **Equivalent weight flow W√θ/δ** | — | **0.581 lb/s** |
| Total-pressure ratio P₆/P₁ | 1.90 | 1.90 |
| Overall efficiency η₁₋₆ | 0.80 | 0.80 |
| Rotative speed N | 36 000 rpm | 52 200 rpm |
| **Equivalent speed N/√θ** | — | 51 176 rpm |
| Impeller tip speed U_t3 | 667 ft/s | 968 ft/s |
| **Equivalent tip speed U_t3/√θ** | — | **949 ft/s** |
| Impeller exit dia. / blades | 4.25 in (10.8 cm) / 15 (30° backsweep), 17 diffuser vanes | — |

The **same physical aerodynamic point** appears at 36 000 rpm in He–Xe and 52 200 rpm in argon
(ratio 1.45 ≈ √(MW_HeXe/MW_Ar)) — the corrected-speed (Mach) similarity that the whole
gas-transfer rests on (§4).

## 2. Corrected-parameter definitions (NASA standard-day referral)

Reference: **T_std = 518.67 °R, P_std = 14.696 psia.**

```
θ = T₁/T_std    δ = P₁/P_std
W_corr (equivalent flow) = W·√θ/δ      N_corr = N/√θ      U_corr = U/√θ
```

Verification at the argon design point (θ=1.0411, δ=1.3779):
W·√θ/δ = 0.785·1.0203/1.3779 = **0.5812 lb/s** ✓ (Table II 0.581); N/√θ = **51 161 rpm** ✓ (0.581…);
U/√θ = **948.7 ft/s** ✓ (949). So the report's "equivalent" *is* standard-day corrected — digitized
x-axes map straight onto NPSS `WcMap`, speed labels onto `NcMap`.

## 3. Map construction

NPSS β-line (R-line) map: `NcMap = speed%/100`, `WcMap = We/We_des` (We_des = 0.581 lb/s, a Table-II
constant), `betaMap` 0 = surge … 1 = max flow (β linear in equivalent flow per speed line). Outputs
`PRmap = P6/P1` and `effMap = η₁₋₆`.

Each speed line's `PR(We)` and `eff(We)` are fit with a **shape-preserving PCHIP** (piecewise-cubic
Hermite), then sampled on the β grid. PCHIP is chosen over a quartic/smoothing spline because those
*overshoot* above the data on the gentle (surge) side of the asymmetric efficiency island and dip at
the flat PR tops — distorting the physical shape. PCHIP is smooth (C¹), monotone-respecting, and
introduces no spurious extrema. Because compressor PR is non-monotonic (rises to a peak near surge,
then falls) and η is single-peaked, there is no single closed form — PCHIP is the right tool (the
turbine, by contrast, has a clean swallowing-law form).

**Digitization (the foundation).** Both PR (Fig 8) and efficiency (Fig 10) come from the user's
manual marker-center WebPlotDigitizer reads (`data/manual_digitize/`), ingested by
`scripts/build/ingest_manual_digitize.py`. This replaced an earlier pixel read that carried a
systematic ~16 % flow-axis scale error; the manual reads were validated three independent ways on
the 100 % line (agreeing to ~0.002 lb/s) and against the printed design dot. Design dot now reads
We = 0.577 (spec 0.581).

**Where the design point lands.** With We_des = 0.581, the design point sits at `NcMap = 1.000`,
`betaMap ≈ 0.606`. Querying the map at (NcMap=1.0, WcMap=1.0) returns **PR = 1.912, η = 0.798** vs
the spec 1.90 / 0.795 — a ~1 % residual because the Table-II design flow (0.581) sits just right of
the plotted design marker (0.577) on the 100 % line. Peak η is ~0.819 at design speed, rising to
~0.83 at 50 % speed (low-power off-design). Normalized `WcMap` range 0.25 … 1.16.

## 4. Gas transfer (argon rig → He-Xe / Kr) — validated

NASA "equivalent" flow/speed are **gas-referenced** (scale as √MW), so they differ between fluids for
the *same* aerodynamic point. But the **normalized** map (`NcMap = N/N_des`, `WcMap = We/We_des`)
cancels the √MW factor — **the normalized map is gas-independent**; the fluid enters only through the
per-gas design scalars used to set the NPSS map scalers (`s_Wc`, `s_Nc`). All candidate fluids are
monatomic (γ = 5/3), so only Reynolds differs.

Validated against the real He-Xe loop (NASA TM X-67989, Fig 9): He-Xe design equivalent flow 0.85 vs
argon 0.581 → ratio 1.46 = √(83.8/39.95) = 1.448 ✓; design speed 36 000 vs 52 200 rpm → 1.45 ✓;
design PR 1.90 = 1.90 ✓. (The design He-Xe MW 83.8 equals krypton's MW, so Kr and He-Xe share scalars.)

## 5. Corrections for real-loop representation

The rig map is a clean, near-design-Reynolds aerodynamic map. Variants layer the real-loop effects.

**Reynolds (NASA TN D-6640).** `Re_U = ρ₁·U_t3·D_t3/μ₁`, design 3.1×10⁶. Loss law
`(1−η)/(1−η_ref) = [(Re_U)_ref/Re_U]ⁿ` with piecewise n (0.06 / 0.09 / 0.20 from high to low Re).
Magnitude η_max 0.813 (design Re) → 0.805 (Re 1.5e6) → 0.792 (7.0e5) → 0.773 (3.4e5), matching the
report's "1.5 pts to 30 %, +2.5 pts to 10 %." `scripts/build/reynolds_correction.py` +
`generate_map_variants.py` → `maps/compressor/variants/{rig_design_Re100, loop_Re48, loop_Re23,
loop_Re11}`. Surge/max-flow shift to lower flow at low Re is *not* applied (no published number;
reserved for loop calibration).

**Build / installation variation (TM X-52826, three BRUs).** Hardware build variation is
**negligible** — He-Xe and Kr plot on common curves and the three units are indistinguishable. The
74 % (BRU rig) vs 80 % (gas loop) compressor spread is an **operating-point effect** (the rig runs
open-throttle at high flow, down the efficiency ridge), already captured by the map — not a hardware
band. Residual hardware+instrument band is small: **±1.5 η pts, ±2 % PR**.

**Heat-soak "indicated" efficiency.** The insulated hot loop soaks heat into the compressor → high
measured discharge temperature → indicated (temperature-rise) efficiency reads low (loop 69–75 % vs
rig 80–82 %). Modeled as a parasitic dT_soak/T₁ ≈ 0.033.

**Calibrated set** (`scripts/build/calibration.py` + `generate_calibrated_variants.py` →
`maps/compressor/calibrated/`):

| variant | peak η (100%) | use |
|---|---|---|
| aero_nominal | 0.821 | true aerodynamic (Reynolds-corrected); model heat-soak in the cycle |
| aero_optimistic / pessimistic | 0.836 / 0.806 | hardware + instrument band |
| indicated_nominal | 0.756 | loop as-instrumented (heat-soak baked in) |
| indicated_pessimistic | 0.741 | loop worst case |

## 6. Validation against hot closed-loop data (NASA TM X-2350)

The BRU compressor was run in the actual hot loop vs the Compressor Research Package (CRP) map:
- CRP map → **80 % at BRU design**; our aero map peak 0.821, design 0.798 ≈ 80 % ✓.
- In-loop BRU compressor → **74.5 %**, decomposed (by the report) into exactly the effects we
  identified: axial clearance (7.5 % vs 4 % → −3.5 % flow), seal leakage (−0.3 pt), and off-design
  open-throttle operation. Rotor power balance closes to **0.67 %** of ideal turbine power.
- Turbine→compressor heat transfer ≈ 1.5 η pts here — small (well-insulated), and a *separate,
  condition-dependent* term that belongs at the cycle level, not baked into the map (the project's
  "Option B": clean aero maps, heat leak modeled in the cycle).

> Reynolds caveat: TM X-2350's compressor Reynolds (design ≈ 251 000) uses a **different length
> scale** than the tip-based Re_U ≈ 3.1×10⁶ in TN D-6640. Same physics, different definition — do
> not conflate the two numbers.

## 7. Integration (verified, runnable)

`scripts/test/test_map_integration.jl` loads the map into the sibling GasCycle solver:
```
read_npss_map(compressor_argon.map) → to_performance_map → query(Nc=1.0, Wc=1.0)
  => PR = 1.912, eta = 0.798     Nc_axis [0.5 … 1.0]   Wc_axis 0.253 … 1.156
```
The β-line NPSS `Table` map is consumable with no solver code changes; `scale_map` then anchors it
to the chosen design point for any working fluid.
