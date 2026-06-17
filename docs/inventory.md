# BRU Document Inventory & Data Triage

Corpus: 23 NASA reports (1967–1972) on the Brayton Rotating Unit (BRU) and its
4.25-in centrifugal compressor and 4.97-in radial-inflow turbine. Goal: reconstruct
compressor & turbine performance maps for an NPSS (or similar) cycle model of the
recuperated closed Brayton cycle, with emphasis on off-design.

Roles: **[ANCHOR]** geometry/design definition · **[MAP]** primary map source (rig data)
· **[CORR]** correction data (Reynolds, clearance) · **[GAS]** working-fluid variant for
validation · **[SYS]** system-level validation · **[REF]** background.

## Compressor (4.25-in / 10.8-cm sweptback-bladed centrifugal)

| Report | Role | What it gives |
|---|---|---|
| 1970 Overall Performance in Argon of 4.25-Inch…Compressor | **[MAP] PRIMARY** | Full map in argon. Fig 8 = PR vs W_eq (6 speed lines 50–100% + surge line + design pt). Fig 9 = ΔT/T₁ vs W_eq. Fig 10 = η vs W_eq. Fig 11 = work factor. Tables I (HeXe design pt) & II (argon design pt + equivalent values). Inlet held at 20 psia/540°R. |
| 1972 Reynolds Number Effect on Overall Performance of 10.8-cm…Compressor | **[CORR]** | Reynolds-number correction for the compressor map (low-density / off-design). |
| 1971 Preliminary Performance of the Brayton 4.25-inch Radial Compressor…Helium Xenon | **[GAS]** | Same compressor in HeXe (MW≈83.8) over TIT 1200–1600°F. Validates argon→HeXe corrected-coordinate transfer. Also has argon & krypton refs. |
| 1967 Compressor Research Package Final Report | **[ANCHOR]** | Original compressor aerodynamic design (geometry, design velocity diagrams, design intent). |

## Turbine (4.97-in radial-inflow)

| Report | Role | What it gives |
|---|---|---|
| 1969 Cold Performance Evaluation of 4.97-Inch Radial-Inflow Turbine… | **[MAP] PRIMARY** | Full cold-flow map in argon. 0–110% design speed, PR 1.33–2.01. Mass flow, torque, efficiency, specific work; Reynolds sweep (Re 34,950–175,800; design Re 76,200). |
| 1969 Effect of Axial Running Clearance on Performance of Two Brayton Cycle Radial Inflow Turbines | **[CORR]** | Tip/axial-clearance correction for turbine efficiency & flow. |
| 1971 Preliminary Performance of a 4.97-Inch Radial Turbine…Helium-Xenon | **[GAS]** | Turbine in HeXe (MW≈83.5–83.8). Validates argon→HeXe transfer. (Two near-identical copies in corpus.) |
| 1972 Recent Radial Turbine Research at the NASA Lewis Research Center | **[REF]** | Survey/context for the radial turbine program; possible supplemental design data. |
| 1968 Turbine Research Package Final Report | **[ANCHOR]** | Original turbine aerodynamic design (geometry, design point, velocity diagrams). |

## System-level / BRU integration

| Report | Role | What it gives |
|---|---|---|
| 1972 Design and Fabrication of the Brayton Rotating Unit (391 p) | **[ANCHOR]** | Definitive BRU mechanical/aero design & geometry. Large scanned report. |
| 1969 The Design and Fabrication of the BRU on Roller Bearings (232 p) | **[ANCHOR]** | Earlier BRU build; bearing/rotor detail. Large scan. |
| 1971 Turbine and Compressor Performance of a BRU During Hot Closed-Loop Operation | **[SYS]** | Both machines matched, hot closed-loop — direct end-to-end map check. |
| 1970 Experimental Performance characteristics of three identical brayton rotating units | **[SYS]** | Unit-to-unit scatter → uncertainty band on the maps. |
| 1970 Experimental Performance of a 2-15 kW Brayton Power System…Krypton | **[SYS/GAS]** | Krypton system operating points. |
| 1970 Experimental Performance of a 2-15 kW Brayton Power System…Helium and Xenon | **[SYS/GAS]** | HeXe system operating points. |
| 1971 Performance of a Brayton-Cycle Power Conversion System…Helium-Xenon (×2 copies) | **[SYS/GAS]** | HeXe power-conversion-system performance. |
| 1970 Preliminary Performance of a Brayton-Cycle Gas Loop…Krypton (TIT 1200–1600°F) | **[SYS/GAS]** | Krypton gas-loop performance vs TIT. |
| 1970 Performance of the electrically-heated 2 to 15 kWe Brayton power system | **[SYS]** | Electrically-heated system points. |
| 1970 Mechanical Performance of a 2- to 10-kW Brayton Rotating Unit | **[REF]** | Mechanical (bearings, windage) — parasitics for cycle model. |
| 1970 Steady-state analysis of a brayton space power system | **[REF]** | NASA's own cycle model — methodology cross-check & expected operating points. |
| 1972 Predictability of Brayton Electrical Power System Performance | **[REF]** | How well component maps predicted system performance — validation target. |

## Key cross-cutting facts established so far
- All working fluids (argon, krypton, He–Xe) are **monatomic, γ = 5/3**. Argon was used as
  a rig surrogate specifically because it matches the design HeXe specific-heat ratio.
  ⇒ corrected-coordinate maps transfer across fluids with only a **Reynolds-number** correction.
- Compressor "equivalent" quantities use **NASA standard day (518.67 °R, 14.696 psia)** referral.
  Verified: argon design W=0.785 lb/s @ 20.25 psia/540°R → W√θ/δ = 0.581 lb/s (matches Table II).
- Compressor map speed lines: U_tip,eq = 949/854/759/664/569/475 ft/s = 100/90/80/70/60/50% design.
