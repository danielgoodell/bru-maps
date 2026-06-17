"""
Generate Reynolds-corrected compressor map VARIANTS for different operating situations.

Takes the base argon rig map (maps/compressor/compressor_map_gridded.csv, which is at reference
Re_U ~ 3.1e6) and produces one corrected map per scenario, so the cycle model can use the
map that matches the real test-loop density/Reynolds regime.

Each variant -> a gridded CSV and an NPSS .map file in maps/compressor/variants/, plus one overlay
plot comparing efficiency (the dominant Reynolds effect) across variants.

Scenarios here are bracketed by Re_U FRACTION of design (3.1e6). The actual loop Re_U for
the He-Xe / Kr power-system tests will be pinned in Phase 5 from the system reports and the
verified gas viscosities; until then these are illustrative-but-physical brackets.
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
BETA = np.linspace(0.0, 1.0, 11)

# name -> Re_U (absolute).  Edit/extend as loop conditions are pinned.
SCENARIOS = {
    "rig_design_Re100": 3.10e6,   # argon rig, near design Re (= base map)
    "loop_Re48":        1.50e6,   # ~48% design Re (higher-pressure loop point)
    "loop_Re23":        0.70e6,   # ~23% design Re (mid loop point)
    "loop_Re11":        0.34e6,   # ~11% design Re (lowest Re tested)
}


def load_base():
    return pd.read_csv(os.path.join(MAPS, "compressor_map_gridded.csv"))


def write_npss(df, Re_U, path, name):
    Nc = sorted(df.NcMap.unique())
    S = float(rc.loss_scale(Re_U))
    eta_max = 1 - rc.LOSS_REF * S
    hdr = (f"// BRU 4.25-in compressor map VARIANT '{name}'  Re_U={Re_U:.2e} "
           f"({100*Re_U/RE_DESIGN:.0f}% design)\n"
           f"// Reynolds-corrected from NASA TM X-2129 rig map via NASA TN D-6640 loss law.\n"
           f"// loss-scale S={S:.4f}  -> max adiabatic eff ~ {eta_max:.3f} "
           f"(rig {1-rc.LOSS_REF:.3f}).  PR knock-down applied; flow/surge shift NOT (calibrate).\n"
           f"// Independents NcMap, betaMap; outputs PRmap, effMap, WcMap (Wc=We/{COMP.We_des_lbps}).\n")

    def table(col, out):
        lines = [f"Table S_map.{out}(real NcMap, real betaMap) {{"]
        for nc in Nc:
            d = df[df.NcMap == nc].sort_values("betaMap")
            betas = ", ".join(f"{b:.2f}" for b in d.betaMap)
            vals = ", ".join(f"{v:.4f}" for v in d[col])
            lines += [f"   NcMap = {nc:.3f} {{",
                      f"      betaMap = {{ {betas} }}",
                      f"      {out} = {{ {vals} }}", "   }"]
        lines.append("}")
        return "\n".join(lines)

    with open(path, "w") as f:
        f.write(hdr + "\n" + "\n\n".join(
            table(c, o) for c, o in [("PRmap", "PRmap"), ("effMap", "effMap"),
                                     ("WcMap", "WcMap")]) + "\n")


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
