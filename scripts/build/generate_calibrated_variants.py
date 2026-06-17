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


def write_npss(df, path, name, note):
    Nc = sorted(df.NcMap.unique())
    hdr = (f"// BRU 4.25-in compressor map - calibrated variant '{name}'\n"
           f"// {note}\n"
           f"// Independents NcMap, betaMap; outputs PRmap, effMap, WcMap (Wc=We/{COMP.We_des_lbps}).\n")

    def table(col, out):
        lines = [f"Table S_map.{out}(real NcMap, real betaMap) {{"]
        for nc in Nc:
            d = df[df.NcMap == nc].sort_values("betaMap")
            betas = ", ".join(f"{b:.2f}" for b in d.betaMap)
            vals = ", ".join(f"{v:.4f}" for v in d[col])
            lines += [f"   NcMap = {nc:.3f} {{", f"      betaMap = {{ {betas} }}",
                      f"      {out} = {{ {vals} }}", "   }"]
        lines.append("}")
        return "\n".join(lines)

    with open(path, "w") as f:
        f.write(hdr + "\n" + "\n\n".join(table(c, c) for c in ("PRmap", "effMap", "WcMap")) + "\n")


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
