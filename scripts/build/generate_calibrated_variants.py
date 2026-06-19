"""
Generate the calibrated compressor map SET for real-loop cycle modeling.

From the base argon rig map (maps/compressor/compressor_map_gridded.csv) produce:
  aero_nominal        - Reynolds-corrected rig (true aerodynamic) efficiency  [Option B base]
  aero_optimistic     - aero + hardware/instrument band (+1.5 pts eff, +2% PR)
  aero_pessimistic    - aero - band
  indicated_nominal   - heat-soak 'indicated' efficiency (peak ~0.73-0.75)     [Option A]
  indicated_pessimistic - indicated - band (loop worst case)

Each -> maps/compressor/calibrated/<name>.map (NPSS) + .csv, plus one comparison plot.

Evidence base: NASA TM X-52826 (three identical BRUs) shows hardware/gas variation is
negligible and the installation efficiency spread is an operating-point effect already in the
map; NASA TM X-67989 gives the 69-75% indicated band and the heat-soak cause. See
docs/compressor.md (§5 corrections, §6 validation).
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from corrected_params import COMP
import calibration as cal

HERE = os.path.dirname(__file__)
MAPS = os.path.join(HERE, "..", "..", "maps", "compressor")
OUT = os.path.join(MAPS, "calibrated")
os.makedirs(OUT, exist_ok=True)

ALPHA_LEVELS = [0.0, 1.0]   # two identical sub-maps + linear alpha interp (fixed-geometry BRU)


def _rline_des(df):
    """Design R-line: where the 100% line crosses WcMap = 1.0 (unchanged from the base map)."""
    d = df[df.NcMap == 1.0].sort_values("WcMap")
    return float(np.interp(1.0, d.WcMap, d.RlineMap))


def write_npss(df, path, name, note):
    Nc = sorted(df.NcMap.unique())
    hdr = (f"// BRU 4.25-in compressor map - calibrated variant '{name}'\n"
           f"// {note}\n"
           f"// NPSS R-line map ('Subelement CompressorRlineMap S_map'). Independents alphaMap "
           f"(fixed-geom: 0.0), NcMap, RlineMap (1.0=surge..3.0=max flow); tables TB_Wc, TB_PR, TB_eff "
           f"(leaves WcMap=We/{COMP.We_des_lbps}, PRmap, effMap).\n"
           f"// Surge/stall line = RlineMap 1.00 column (RlineStall below); efficiency-band "
           f"calibration does not move the surge flow.\n")

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
    base = pd.read_csv(os.path.join(MAPS, "compressor_map_gridded.csv"))
    eff0, PR0 = base.effMap.values, base.PRmap.values
    ind0 = cal.indicated_efficiency(eff0, PR0)

    variants = {}
    variants["aero_nominal"] = (eff0, PR0, "Reynolds-corrected rig aero efficiency (true aero).")
    e, p = cal.apply_band(eff0, PR0, +1)
    variants["aero_optimistic"] = (e, p, "aero + hardware/instrument band (+1.5pts eff, +2% PR).")
    e, p = cal.apply_band(eff0, PR0, -1)
    variants["aero_pessimistic"] = (e, p, "aero - hardware/instrument band.")
    variants["indicated_nominal"] = (ind0, PR0,
        "heat-soak 'indicated' efficiency (insulated hot loop); peak ~0.73-0.75 matches measured.")
    e, p = cal.apply_band(ind0, PR0, -1)
    variants["indicated_pessimistic"] = (e, p, "indicated - band (loop worst case).")

    for name, (eff, PR, note) in variants.items():
        df = base.copy()
        df["effMap"] = np.clip(eff, 0.05, 0.95)
        df["PRmap"] = PR
        df.to_csv(os.path.join(OUT, f"{name}.csv"), index=False)
        write_npss(df, os.path.join(OUT, f"{name}.map"), name, note)

    # comparison plot: 100% speed-line efficiency for each variant
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    m = base.NcMap == 1.0
    We = base[m].WcMap * COMP.We_des_lbps
    styles = {"aero_nominal": ("C0", "-"), "aero_optimistic": ("C0", ":"),
              "aero_pessimistic": ("C0", "--"), "indicated_nominal": ("C3", "-"),
              "indicated_pessimistic": ("C3", "--")}
    for name, (eff, PR, _) in variants.items():
        c, ls = styles[name]
        ax.plot(We, np.clip(eff, 0.05, 0.95)[m.values], ls, color=c, lw=1.8, label=name)
    # reference markers from the three-BRU data
    ax.axhline(0.80, color="k", lw=0.6, ls=":"); ax.text(0.46, 0.806, "gas-loop / ref ~0.80", fontsize=7)
    ax.axhspan(0.69, 0.75, color="orange", alpha=0.12)
    ax.text(0.46, 0.70, "measured indicated band 69-75%", fontsize=7, color="darkorange")
    ax.set(xlabel="Equiv. weight flow W√θ/δ (argon), lb/s", ylabel="efficiency",
           title="Compressor 100% line — calibrated variant set (aero vs indicated + band)")
    ax.grid(alpha=0.3); ax.legend(fontsize=8, loc="lower left")
    fig.tight_layout(); fig.savefig(os.path.join(MAPS, "compressor_calibrated_variants.png"), dpi=120)

    print(f"Wrote {len(variants)} calibrated maps to maps/compressor/calibrated/ + comparison plot.")
    for name, (eff, PR, _) in variants.items():
        pk = float(np.nanmax(np.clip(eff, 0.05, 0.95)[m.values]))
        print(f"  {name:24s} peak-100%-eff = {pk:.3f}")


if __name__ == "__main__":
    main()
