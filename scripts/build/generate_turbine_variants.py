"""
Generate turbine map variants for different operating situations from the base cold-rig map
(maps/turbine/turbine_map_gridded.csv = aero, design Re, ~2% clearance).

  aero_nominal       - cold-rig aerodynamic map (design Re, reference clearance)
  loop_clearance12   - in-loop axial clearance 12% (TM X-2350), design Re
  loop_clr12_lowRe   - 12% clearance + low-power Reynolds (Re ~46 700, low end of loop)

Each -> maps/turbine/variants/<name>.map + .csv, plus a comparison plot.
Corrections from turbine_corrections.py (NASA TN D-5090 Reynolds, TM X-52552 clearance).
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import turbine_corrections as tc

HERE = os.path.dirname(__file__)
MAPS = os.path.join(HERE, "..", "..", "maps", "turbine")
OUT = os.path.join(MAPS, "variants")
os.makedirs(OUT, exist_ok=True)
Wc_DES = 0.486


def write_npss(df, path, name, note):
    Nc = sorted(df.NpMap.unique())
    hdr = (f"// BRU 4.97-in turbine map variant '{name}'\n// {note}\n"
           f"// NPSS turbine map ('Subelement TurbinePRmap S_map'): independents NpMap, PRmap "
           f"(total-to-static p1'/p3); tables TB_Wp, TB_eff output WpMap, effMap.\n")

    def table(tb_name, leaf):
        lines = [f"   Table {tb_name}(real NpMap, real PRmap) {{"]
        for nc in Nc:
            d = df[df.NpMap == nc].sort_values("PRmap")
            prs = ", ".join(f"{p:.3f}" for p in d.PRmap)
            vals = ", ".join(f"{v:.4f}" for v in d[leaf])
            lines += [f"      NpMap = {nc:.3f} {{", f"         PRmap = {{ {prs} }}",
                      f"         {leaf} = {{ {vals} }}", "      }"]
        lines += ['      NpMap.interp = "lagrange3";  NpMap.extrap = "linear";',
                  '      PRmap.interp = "lagrange3";  PRmap.extrap = "linear";',
                  "      extrapIsError = 0;",
                  "      printExtrap   = 0;"]
        lines.append("   }")
        return "\n".join(lines)

    sub = [
        "Subelement TurbinePRmap S_map {",
        "   PRmapDes = 1.690;",
        "   NpMapDes = 1.000;",
        "",
        table("TB_Wp", "WpMap"),
        "",
        table("TB_eff", "effMap"),
        "}",
    ]
    with open(path, "w") as f:
        f.write(hdr + "\n" + "\n".join(sub) + "\n")


def main():
    base = pd.read_csv(os.path.join(MAPS, "turbine_map_gridded.csv"))
    variants = {}
    variants["aero_nominal"] = (base.effMap.values, base.WpMap.values,
                                "cold-rig aero, design Re, ~2% clearance")
    e, w = tc.correct_clearance(base.effMap.values, base.WpMap.values, 12.0)
    variants["loop_clearance12"] = (e, w, "in-loop 12% axial clearance (TM X-2350), design Re")
    e2 = tc.correct_efficiency_reynolds(e, 46700.0)
    variants["loop_clr12_lowRe"] = (e2, w, "12% clearance + low-power Re~46700")

    for name, (eff, wc, note) in variants.items():
        df = base.copy()
        df["effMap"] = np.clip(eff, 0.05, 0.97)
        df["WpMap"] = wc
        df.to_csv(os.path.join(OUT, f"{name}.csv"), index=False)
        write_npss(df, os.path.join(OUT, f"{name}.map"), name, note)

    fig, ax = plt.subplots(figsize=(8, 5.5))
    m = base.NpMap == 1.0
    PR = base[m].PRmap
    for name, (eff, wc, _) in variants.items():
        ax.plot(PR, np.clip(eff, 0.05, 0.97)[m.values], "-o", ms=3, label=name)
    ax.axhline(0.913, color="k", lw=0.6, ls=":"); ax.text(1.45, 0.916, "cold-rig design η_t 0.913", fontsize=7)
    ax.axhline(0.883, color="r", lw=0.5, ls=":"); ax.text(1.45, 0.886, "in-loop η(1→3)≈0.88 @12% clr", fontsize=7, color="r")
    ax.set(xlabel="(p1'/p3) total-to-static PR", ylabel="total efficiency η_t (1→3)",
           title="Turbine 100% line — clearance / Reynolds variants")
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(MAPS, "turbine_variants.png"), dpi=120)

    print(f"Wrote {len(variants)} turbine variant maps to maps/turbine/variants/ + plot.")
    for name, (eff, wc, _) in variants.items():
        print(f"  {name:20s} peak η = {float(np.nanmax(np.clip(eff,0.05,0.97)[m.values])):.3f}")


if __name__ == "__main__":
    main()
