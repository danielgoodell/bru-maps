# Turbine map — 4.97-in radial-inflow turbine

Primary source: **NASA TN D-5090 (1969)**, *Cold Performance Evaluation of the 4.97-Inch
Radial-Inflow Turbine* (argon rig; Fig 8 flow-PR, Fig 11a static-η and 11b total-η vs blade-jet
speed ratio). The turbine analog of the compressor TM X-2129.

Deliverables: `maps/turbine/turbine_argon.map` (total-to-static PR axis) and
`maps/turbine/turbine_argon_tt.map` (total-to-total PR axis) + gridded CSVs, validation PNG, and the
variant set. Built by `scripts/build/build_turbine_map.py` from `data/digitized/turbine_*_argon.csv`.

---

## 1. Design point (Table I)

| Parameter | Design | Experimental |
|---|---|---|
| Total efficiency η_t (1→3, rotor exit) | 0.897 | 0.913 |
| Static efficiency η_s (1→3) | 0.850 | 0.872 |
| Overall total efficiency (1→4, diffuser exit) | 0.884 | 0.894 |
| Equivalent mass flow εW√θcr/δ, lb/s | 0.4860 | 0.485 |
| Blade-jet speed ratio ν at design | 0.70 | 0.70 |
| Design Reynolds number | 76 200 | — |
| Tip diameter | 4.97 in (12.62 cm) | — |

Design total-to-static pressure ratio (p₁'/p₃)_eq ≈ **1.69** (Fig 8 design marker).

## 2. Equivalent (corrected) parameters — critical-velocity form

Standard 518.67 °R / 14.70 psia; subscript cr = Mach-1 (critical). `δ = p₁'/p*`, `θcr` = inlet
critical-velocity-squared ratio, `ε` = γ-dependent gas-referral factor (to air).
```
equivalent mass flow = ε·W·√θcr/δ      equivalent speed = N/√θcr
blade-jet speed ratio ν = U_t/V_j       (rotor-inlet tip speed / ideal jet speed)
```
For monatomic gases (Ar, Kr, He-Xe; γ = 5/3) ε is common, so the equivalent map is the same for all
three fluids — only Reynolds differs (same conclusion as the compressor; §6 validates it for He-Xe).

## 3. Map construction

NPSS β-line map: `NcMap = speed%/100`, `WcMap = Wc/Wc_des` (Wc_des = 0.486 lb/s, Table-I constant),
`betaMap` 0 = low PR … 1 = high PR. Two ideas:

**Flow (Fig 8).** Each speed line is fit independently to a three-parameter swallowing law
`Wc(PR) = W_ch·(1 − PR^−a)^b` — monotonic increasing, Wc(PR=1)=0, saturates toward W_ch. The free
exponent `b` is essential: every 2-parameter form (ellipse `b=0.5`, exponential, Michaelis) leaves a
systematic `+ − − − + +` residual *arch* on every line — the measured curves bend harder than any
fixed-curvature form. The free `b` removes that and drops residuals to the noise floor (per line
0.2–1.2×10⁻³ lb/s, ~0.1–0.2 %; 110 % line highest at the high-PR end). Fit per-line, not as a
smooth-in-speed family — the lines differ enough in curvature that a global model reintroduces
structured residuals; independent fits each reach the floor and stay ordered (the data don't cross).
The map never extrapolates: each line uses its own measured PR range, off-design interpolates the
6-line (30–110 %) grid.

**Efficiency (Fig 11).** Radial-turbine efficiency collapses onto a single curve of η vs ν. The
user digitized **both** the static (Fig 11a) and total (Fig 11b) curves per speed line; the build
**pools every speed line's markers** into one (ν, η) cloud and fits one polynomial each — the
measured collapse itself. Efficiency at any (Nc%, PR) follows from
```
ν(Nc%, PR) = 0.70·(Nc%/100)·√(1 − 1.69^−0.4) / √(1 − PR^−0.4)      (γ = 5/3, exp 0.4)
```
then η = η_total(ν). At design: ν(100 %, 1.69) = 0.700, η_total(0.70) = 0.916 (spec 0.913, within
read precision).

**Design flow now coincides with the 100 % line.** In the older pixel read the Fig-8 design dot sat
well above the faired 100 % line (a +2.5–3.3 % figure-vs-table offset that forced an awkward anchor
and a kink). The user's manual reads remove it: the fitted 100 % line at PR 1.69 gives Wc = 0.485 vs
the Table-I 0.486 — **−0.2 %**. So design flow and design PR coincide on the map: querying at
(NcMap=1.0, WcMap=1.0) returns **PR = 1.694, η = 0.916** (vs design 1.69 / 0.913). No anchor, no kink.

## 4. PR basis: total-to-static ↔ total-to-total (and the efficiency-definition trap)

The report's native PR axis is **inlet-total / exit-static** (p₁'/p₃) — it sets V_j and the flow. An
NPSS/GasCycle turbine element uses **total-to-total** expansion ratio Pt_in/Pt_out. A turbine's
actual work is one physical quantity, normalized two ways:
```
η_tt = Dh_act/Dh_ideal(→ exit TOTAL p)      η_ts = Dh_act/Dh_ideal(→ exit STATIC p)
```
The total-to-static ideal work is larger (by the exit kinetic energy), so **η_ts < η_tt always**.
Design: η_tt = 0.916, η_ts = 0.873 (gap 0.043 = exit-KE share). Correspondingly PR_ts > PR_tt;
the exit-KE is recovered directly from the two digitized ν-curves (no assumed exit Mach):
```
PR_tt = [ 1 − (η_ts/η_tt)·(1 − PR_ts^−0.4) ]^(−1/0.4)        design: PR_ts 1.69 → PR_tt 1.645
```
`build_turbine_map.py` emits both `turbine_argon.map` (PR_ts axis) and **`turbine_argon_tt.map`**
(PR_tt axis — the one to pair with a total-to-total cycle element), plus `turbine_map_tt_gridded.csv`
(every point with ν, PR_ts, PR_tt, η_tt, η_ts).

**The trap (why turbine exit temperature comes out wrong).** Exit total temperature is
`T3' = T1'·(1 − η·(1 − PR^−0.4))`, correct **only if η and PR share a basis**. Pair η_tt with PR_tt,
OR η_ts with PR_ts — both give the same exit temp. *Mixing* them mis-predicts T3'. A classic symptom:
a deck on a total-to-static PR axis that needs its efficiency "tuned down" (e.g. 0.87 → 0.81) to fix
exit temperature is really mismatching bases — it should use η_ts (~0.84 in-loop), not a
total-to-total 0.87. Using the **tt-map** removes the guesswork: it returns the operating-point
efficiency on a stated basis, so exit temperature is right by construction.

## 5. Corrections for real-loop representation

**Reynolds (TN D-5090 Fig 14).** `Re = W/(μ·r_t)`, design 76 200. Total-efficiency loss law
`(1−η)/(1−η_ref) = (Re_ref/Re)ⁿ` with **piecewise n = 0.30 (Re < 76 200) / 0.18 (Re ≥ 76 200)**,
continuous at design — reproduces the Fig-14 anchors η_t 0.890 / 0.913 / 0.925 at Re
34 950 / 76 200 / 175 800 exactly. (Static eff is less Re-sensitive: n ≈ 0.20 / 0.13.) Fig 14 is a
faired cross-plot with no data points; anchors are text-stated. `scripts/build/turbine_corrections.py`.

**Axial clearance (TM X-52552, 5-in turbine).** η × (1 − 0.00333·Δclr%), flow × (1 − 0.002·Δclr%),
reference ~2 % (cold rig). In-loop 12 % → η(1→3) ×0.967, flow ×0.980.

**Variant set** (`scripts/build/generate_turbine_variants.py` → `maps/turbine/variants/`):
aero_nominal (peak η 0.916), loop_clearance12 (0.885), loop_clr12_lowRe (0.867).

## 6. Cross-validation

**vs in-loop data (TM X-2350).** Research-package overall (1→4) η at BRU equivalent design = 0.894;
in-loop turbine = 0.860. Our clearance model on the overall η: 2 → 12 % gives 0.894 → 0.864 ≈
measured 0.860 ✓ — the bulk of the in-loop deficit is axial clearance (the BRU ran much looser than
the cold rig). The report also warns the in-loop "indicated" efficiency is corrupted by high-temp
heat loss (adiabatic assumption) — another reason to keep maps aerodynamic and handle heat at the
cycle level.

**Gas transfer argon → He-Xe (TM X-67996).** Three checks pass: (1) efficiency physics is
gas-independent — the torque parameter (η/ν) collapses to one curve; our argon η/ν overlays the
He-Xe points to RMS ≈ 0.012, design η/ν = 1.308 vs He-Xe reported 1.305; (2) the He-Xe design point
computes to ν = 0.694 ≈ 0.70 from first principles (PR_ts 1.785, TIT 2060 °R, 36 000 rpm); (3) the
flow function matches at an implied He-Xe turbine-inlet pressure of 41.6 psia (expected ~41–42).
He-Xe loop η ~91 % (incl. diffuser) vs argon 89.4 % overall total-to-total, consistent.

> **He-Xe usage caveat.** The efficiency *physics* (η vs ν) is universal, but the map's
> (NcMap, PR) → ν formula bakes in argon's design PR (1.69). He-Xe's design PR (~1.785) differs
> because PR ↔ ν depends on cp·T₁ (gas *and* temperature). Querying the argon η-grid directly at the
> He-Xe design gives ν = 0.670 → η 0.912 instead of the true ν = 0.694 → η 0.914 (−0.3 eff-pt at the
> flat peak; grows modestly off-design). For best He-Xe accuracy, compute ν from actual conditions
> and look up the universal η(ν), or regenerate the η-grid with PR_des ≈ 1.785. The compressor has
> no analog — its speed/flow normalization absorbs the gas referencing.

## 7. Integration (verified, runnable)

`scripts/test/test_turbine_integration.jl`:
```
read_npss_map(turbine_argon.map) → to_performance_map → query(Nc=1.0, Wc=1.0)
  => PR = 1.694, eta = 0.916     Nc_axis [0.3 … 1.1]
```
Loads into the GasCycle Turbine off-design element with no code changes. For a total-to-total cycle
element, load `turbine_argon_tt.map` instead (design query PR_tt ≈ 1.645).
