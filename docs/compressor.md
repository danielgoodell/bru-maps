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

NPSS **R-line map** — a native `Subelement CompressorRlineMap S_map { … }` block holding three tables
`TB_Wc`, `TB_PR`, `TB_eff` whose leaf arrays are `WcMap = We/We_des` (We_des = 0.581 lb/s, a Table-II
constant), `PRmap = P6/P1`, and `effMap = η₁₋₆`. Each is a function of `(alphaMap, NcMap, RlineMap)`:
`NcMap = speed%/100`; `RlineMap` = compressor R-line running **1.0 (surge) … 3.0 (max flow)** (linear in
equivalent flow per speed line); `alphaMap` is the variable-geometry angle. The BRU is fixed-geometry,
so the map is alpha-independent — but NPSS still interpolates the alpha axis, so each table carries
**two identical sub-maps at `alphaMap = 0` and `1` with linear alpha interpolation**: whatever alpha the
element passes returns the same map, with no need to pin `alphaMap` to 0 (a single sub-map would require
the element to land exactly on it). The design scalars `alphaMapDes`, `NcMapDes`, `RlineMapDes`,
`RlineStall` (and the commented `ReDes` for the Reynolds socket) are member declarations inside the
`S_map` block. Each table then sets its **own** per-axis interpolation/extrapolation methods, declared
*inside the table* after its alpha sub-maps (`lagrange3` cubic in `NcMap`/`RlineMap`, `linear` in
`alphaMap`; `linear` extrapolation, `extrapIsError = 0`, `printExtrap = 0`) — matching the native NPSS
`.map` layout.

**Surge line.** In an R-line map the surge line is not a separate curve — it is the lowest R-line
(lowest corrected flow at each speed). The `S_map` block declares `RlineStall = 1.0`, which is how NPSS
locates the stall boundary, holds the operating `RlineMap ≥ RlineStall`, and computes the stall margins
**SMW** (constant flow) and **SMN** (constant speed). pyCycle's `CompressorRlineMap` uses the identical
convention (`RlineStall = 1.0`, minimum `RlineMap = 1.0`). The build cross-checks
the `RlineMap = 1.0` column against the independently digitized Fig-8 surge line
(`comp_surge_argon.csv`): the two agree to RMS ≈ 0.036 / max ≈ 0.062 in PR, the surge column sitting
just inside the faired line (it is the rig's last stable test point, a hair short of true surge — a
slightly conservative, not optimistic, stall boundary). The 50 % and 70 % lines carry the largest
offset; 80–100 % agree to ≲ 0.014.

Each speed line's `PR(We)` and `eff(We)` are fit with a **shape-preserving PCHIP** (piecewise-cubic
Hermite), then sampled on the R-line grid. PCHIP is chosen over a quartic/smoothing spline because those
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
`RlineMap ≈ 2.21` (set `S_map.RlineMapDes = 2.21`). Querying the map at (NcMap=1.0, WcMap=1.0) returns **PR = 1.912, η = 0.798** vs
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
report's "1.5 pts to 30 %, +2.5 pts to 10 %." Surge/max-flow shift to lower flow at low Re is *not*
applied (no published number; reserved for loop calibration). Source: `reynolds_correction.py`.

**Where the He-Xe operating point sits (`RNI`).** The `RNI = 1` reference is the *argon rig* design
`Re_U` (3.03×10⁶) — **not** the He-Xe mission point. Pinned from the 1971 He-Xe compressor report
(TM X-67846, Table I: inlet 13.5 psia / 540 °R, 36 000 rpm, U_tip 667 ft/s) and the 1970 He-Xe 2–15 kWe
system report (net power = 5.1 kW @ 26 psia → 10.5 kW @ 44.4 psia compressor-discharge, 1600 °F TIT):
- He-Xe **~6 kWe reference design** (13.5 psia inlet): `Re_U ≈ 2.8×10⁶`, **`RNI ≈ 0.92`**.
- He-Xe **~10 kWe net operational point** (≈44 psia discharge → ≈23–24 psia inlet, the *real* design):
  `Re_U ≈ 4.8–5.0×10⁶`, **`RNI ≈ 1.6`** → `s_effRe ≈ 1.006` (a ~+0.5 η-pt **favorable** correction).

So the cold argon rig (`RNI = 1`) was deliberately set at the **low end** of the He-Xe operating
Reynolds range; the 10 kWe point runs ~1.6× higher, where the embedded socket *raises* efficiency
slightly rather than penalizing it. (Inlet pressure ≈ discharge/PR with PR 1.82–1.9; inlet T taken at
the design 540 °R/80 °F — the system swept 60–120 °F; discharge taps may be static, so treat `RNI` as
≈ 1.5–1.7.) Computed via `reynolds_correction.reynolds_U`.

Two ways to consume this, from the same loss law:
- **NPSS-native (S_Re socket), embedded in the base map.** `build_compressor_map.py` embeds a nested
  `Subelement CompressorReynoldsEffects S_Re { … }` (the S_Re socket, type `COMPRESSOR_REYNOLDS_EFFECTS`,
  returns `s_WcRe`, `s_effRe`) directly inside `compressor_argon.map`'s `S_map` block — so the base map
  is self-contained and the correction is a **no-op at the reference Re** (`RNI = 1`). It provides the
  direct multipliers `s_effRe` and `s_WcRe` as functions of `RNI = Re_U/ReDes` (ReDes = 3.03×10⁶). NPSS
  standardizes only the *usage* (`effBase = s_effDes·s_effRe·effMap`, `WcBase = s_WcDes·s_WcRe·WcMap`);
  the subelement internals are implementation-defined, so the `RNI`-indexed lookup-table form is our
  chosen convention. The block is generated by `generate_reynolds_maps.py` (`compressor_s_re_block()`),
  which also writes `compressor_reynolds_correction.csv` as an inspectable table dump. The class name
  follows the NPSS manual's `TurbineReynoldsEffects` example.
  NPSS applies `effBase = s_effDes·s_effRe·effMap` and `WcBase = s_WcDes·s_WcRe·WcMap` — **efficiency
  and flow only; NPSS has no PR Reynolds scalar**, so the small peak-PR knock-down lives only in the
  pre-baked variants. `s_effRe` is anchored at η_ref = 0.813 (`= [1−(1−η_ref)·S]/η_ref`, with S the
  report's deficit law), exact at η_ref and a close uniform-multiplier approximation elsewhere; the
  CSV also carries the deficit S for an S_Re subelement that wants the exact per-point form. `s_WcRe ≡ 1`.
- **Pre-baked variants.** `generate_map_variants.py` →
  `maps/compressor/variants/{rig_design_Re100, loop_Re48, loop_Re23, loop_Re11}` — one full map per
  Reynolds level, for tools that don't apply a run-time correction. These apply the deficit law
  *per point*, so they agree with the `s_effRe` table exactly at η_ref and differ slightly off-peak
  (the uniform-multiplier vs per-point distinction above) — use the variants when you want the exact
  per-point deficit, the S_Re table when you want the native NPSS run-time correction.

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
The R-line NPSS `Table` map is consumable with no solver code changes; `scale_map` then anchors it
to the chosen design point for any working fluid.
