"""
Generate Reynolds-corrected compressor map VARIANTS for different operating situations.

Takes the base argon rig map (maps/compressor/compressor_map_gridded.csv, which is at reference
Re_U ~ 3.1e6) and produces one corrected map per scenario, so the cycle model can use the
map that matches the real test-loop density/Reynolds regime.

Each variant -> a gridded CSV and an NPSS .map file in maps/compressor/variants/, plus one overlay
plot comparing efficiency (the dominant Reynolds effect) across variants.

Scenarios here are bracketed by Re_U FRACTION of design (3.1e6). The He-Xe COMPRESSOR operating
Reynolds is now pinned from the system reports (see docs/compressor.md sec.5): RNI ~0.92 at the
6 kWe reference design and ~1.6 at the 10 kWe net operational point -- i.e. the He-Xe loop runs AT or
ABOVE the rig Re, so these low-Re fractional brackets are illustrative low-power / off-design scenarios,
not the design point. (Kr and the hot turbine Re are still to be pinned.)
Efficiency & PR corrections are quantified from NASA TN D-6640; the flow/surge-to-lower-flow
shift it documents is left for loop-data calibration (not fabricated here).
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from corrected_params import COMP
import reynolds_correction as rc

HERE = os.path.dirname(__file__)
MAPS = os.path.join(HERE, "..", "..", "maps", "compressor")
VAR = os.path.join(MAPS, "variants")
os.makedirs(VAR, exist_ok=True)

RE_DESIGN = 3.1e6
ALPHA_LEVELS = [0.0, 1.0]   # two identical sub-maps + linear alpha interp (fixed-geometry BRU)

# name -> Re_U (absolute).  Edit/extend as loop conditions are pinned.
SCENARIOS = {
    "rig_design_Re100": 3.10e6,   # argon rig, near design Re (= base map)
    "loop_Re48":        1.50e6,   # ~48% design Re (higher-pressure loop point)
    "loop_Re23":        0.70e6,   # ~23% design Re (mid loop point)
    "loop_Re11":        0.34e6,   # ~11% design Re (lowest Re tested)
}


def load_base():
    return pd.read_csv(os.path.join(MAPS, "compressor_map_gridded.csv"))


def _rline_des(df):
    """Design R-line: where the 100% line crosses WcMap = 1.0 (unchanged from the base map)."""
    d = df[df.NcMap == 1.0].sort_values("WcMap")
    return float(np.interp(1.0, d.WcMap, d.RlineMap))


def write_npss(df, Re_U, path, name):
    Nc = sorted(df.NcMap.unique())
    S = float(rc.loss_scale(Re_U))
    eta_max = 1 - rc.LOSS_REF * S
    hdr = (f"// BRU 4.25-in compressor map VARIANT '{name}'  Re_U={Re_U:.2e} "
           f"({100*Re_U/RE_DESIGN:.0f}% design)\n"
           f"// Reynolds-corrected from NASA TM X-2129 rig map via NASA TN D-6640 loss law.\n"
           f"// loss-scale S={S:.4f}  -> max adiabatic eff ~ {eta_max:.3f} "
           f"(rig {1-rc.LOSS_REF:.3f}).  PR knock-down applied; flow/surge shift NOT (calibrate).\n"
           f"// NPSS R-line map ('Subelement CompressorRlineMap S_map'). Independents alphaMap "
           f"(fixed-geom: 0.0), NcMap, RlineMap (1.0=surge..3.0=max flow); tables TB_Wc, TB_PR, TB_eff "
           f"(leaves WcMap=We/{COMP.We_des_lbps}, PRmap, effMap).\n"
           f"// Surge/stall line = RlineMap 1.00 column (RlineStall below). The Reynolds knock-"
           f"down here does NOT shift surge flow, so the stall boundary is unchanged from the base map.\n")

    def table(tb_name, leaf):
        lines = [f"   Table {tb_name}(real alphaMap, real NcMap, real RlineMap) {{"]
        for a in ALPHA_LEVELS:
            lines.append(f"      alphaMap = {a:.3f} {{")
            for nc in Nc:
                d = df[df.NcMap == nc].sort_values("RlineMap")
                rlines = ", ".join(f"{r:.2f}" for r in d.RlineMap)
                vals = ", ".join(f"{v:.4f}" for v in d[leaf])
                lines += [f"         NcMap = {nc:.3f} {{",
                          f"            RlineMap = {{ {rlines} }}",
                          f"            {leaf} = {{ {vals} }}", "         }"]
            lines.append("      }")
        lines += [
            '      alphaMap.interp = "linear";     alphaMap.extrap = "linear";',
            '      NcMap.interp    = "lagrange3";  NcMap.extrap    = "linear";',
            '      RlineMap.interp = "lagrange3";  RlineMap.extrap = "linear";',
            "      extrapIsError = 0;",
            "      printExtrap   = 0;",
        ]
        lines.append("   }")
        return "\n".join(lines)

    sub = [
        "Subelement CompressorRlineMap S_map {",
        "   alphaMapDes = 0.000;",
        "   NcMapDes    = 1.000;",
        f"   RlineMapDes = {_rline_des(df):.3f};",
        "   RlineStall  = 1.000;",
        "   // ReDes    = 3030000.0;   // reference Re_U for the S_Re Reynolds socket (RNI = Re_U/ReDes)",
        "",
        table("TB_Wc", "WcMap"),
        "",
        table("TB_PR", "PRmap"),
        "",
        table("TB_eff", "effMap"),
        "}",
    ]
    with open(path, "w") as f:
        f.write(hdr + "\n" + "\n".join(sub) + "\n")


def main():
    base = load_base()
    fig, ax = plt.subplots(figsize=(8, 5.5))
    speeds = sorted(base.NcMap.unique())
    colors = plt.cm.plasma(np.linspace(0, 0.85, len(SCENARIOS)))

    summary = []
    for k, (name, Re_U) in enumerate(SCENARIOS.items()):
        df = base.copy()
        df["effMap"] = rc.correct_efficiency(base.effMap.values, Re_U)
        df["PRmap"] = rc.correct_PR(base.PRmap.values, Re_U)
        df.to_csv(os.path.join(VAR, f"{name}.csv"), index=False)
        write_npss(df, Re_U, os.path.join(VAR, f"{name}.map"), name)

        # plot the 100% speed-line efficiency for each variant
        d = df[df.NcMap == 1.0]
        We = d.WcMap * COMP.We_des_lbps
        ax.plot(We, d.effMap, "-o", ms=3, color=colors[k],
                label=f"{name}  (Re={Re_U:.2e}, ηmax≈{1-rc.LOSS_REF*float(rc.loss_scale(Re_U)):.3f})")
        summary.append((name, Re_U, float(rc.loss_scale(Re_U)),
                        1 - rc.LOSS_REF * float(rc.loss_scale(Re_U))))

    ax.set(xlabel="Equiv. weight flow W√θ/δ (argon), lb/s",
           ylabel="adiabatic efficiency", title="Compressor 100% speed line — Reynolds variants")
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(MAPS, "compressor_reynolds_variants.png"), dpi=120)

    print(f"{'variant':22s} {'Re_U':>10s} {'loss-scale':>11s} {'eta_max':>8s}")
    for n, re, s, e in summary:
        print(f"{n:22s} {re:10.2e} {s:11.4f} {e:8.3f}")
    print(f"\nWrote {len(SCENARIOS)} variant maps to maps/compressor/variants/ + comparison plot.")


if __name__ == "__main__":
    main()
