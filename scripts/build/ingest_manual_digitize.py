"""
Ingest the user's manual WebPlotDigitizer exports (data/manual_digitize/) into the project's
authoritative digitized inputs (data/digitized/) -- the source of truth for both map builds.

This is STAGE 0 of the pipeline: parse -> validate -> write data/digitized/*.csv (with provenance
headers). Re-run it whenever the manual reads change, then re-run the build scripts.

WPD wide format: row 1 = dataset names (every other column), row 2 = X,Y repeating.
Column/dataset ORDER differs per file (e.g. turbine fig 11b leads with Design), so we key
off row 1 rather than assuming a fixed order. Validation cross-checks each figure's design dot
against the known design point and drops out-of-bounds misclicks; a check figure is written
alongside the exports (data/manual_digitize/manual_digitize_check.png).
"""
import os
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(__file__)
MD = os.path.join(HERE, "..", "..", "data", "manual_digitize")
COMP_DIR = os.path.join(MD, "compressor argon TM X-2129")
TURB_DIR = os.path.join(MD, "turbine argon TN D-5090")
DIG = os.path.join(HERE, "..", "..", "data", "digitized")   # authoritative pipeline inputs

PROV = ("# Source: user MANUAL WebPlotDigitizer reads (data/manual_digitize/), 2026-06-16 -\n"
        "#   geometric marker-center digitization, axis-calibrated to the figure gridlines.\n"
        "#   Supersedes the earlier Claude pixel reads (had a systematic flow-axis scale error);\n"
        "#   regenerate via scripts/ingest_manual_digitize.py. Design dots used as axis sanity checks.\n")


def parse_wpd(path):
    """Return {dataset_name: Nx2 float array} from a WPD wide CSV."""
    with open(path) as f:
        rows = list(csv.reader(f))
    names = rows[0]
    # dataset name sits in the X column; every 2 columns is one dataset
    datasets = {}
    for col in range(0, len(names), 2):
        name = names[col].strip()
        if not name:
            continue
        pts = []
        for r in rows[2:]:
            if col + 1 >= len(r):
                continue
            xs, ys = r[col].strip(), r[col + 1].strip()
            if xs == "" or ys == "":
                continue
            pts.append((float(xs), float(ys)))
        if pts:
            datasets[name] = np.array(pts)
    return datasets


def norm_speed(name):
    """'100%' / '100' -> 100.0 ; non-numeric (e.g. 'Design') -> None."""
    s = name.replace("%", "").strip()
    try:
        return float(s)
    except ValueError:
        return None


def summarize(tag, datasets, xlabel, ylabel, design_spec=None, ybounds=None):
    print(f"\n==== {tag} ====")
    print(f"  axes: X={xlabel}  Y={ylabel}")
    for name, arr in datasets.items():
        x, y = arr[:, 0], arr[:, 1]
        flags = []
        # physical bound checks
        if ybounds is not None:
            lo, hi = ybounds
            if y.min() < lo or y.max() > hi:
                bad = arr[(y < lo) | (y > hi)]
                flags.append(f"OUT-OF-BOUNDS Y: {[tuple(np.round(p,3)) for p in bad]}")
        # monotonic-x ordering check (speed lines should be sorted in flow/PR)
        spd = norm_speed(name)
        if spd is not None and len(x) > 2:
            if not (np.all(np.diff(x) > 0) or np.all(np.diff(x) < 0)):
                order = np.argsort(x)
                flags.append(f"X not monotonic as listed (n={len(x)})")
        n = len(arr)
        print(f"  {name:>8}: n={n:2d}  X[{x.min():.3f},{x.max():.3f}]  Y[{y.min():.3f},{y.max():.3f}]"
              + ("   <-- " + " | ".join(flags) if flags else ""))
    if design_spec is not None:
        for key in ("Design", "design"):
            if key in datasets:
                dp = datasets[key][0]
                dx, dy = design_spec
                print(f"  DESIGN sanity: read=({dp[0]:.3f},{dp[1]:.3f})  spec=({dx:.3f},{dy:.3f})  "
                      f"dX={dp[0]-dx:+.3f} dY={dp[1]-dy:+.3f}")


def write_speedlines(datasets, path, xcol, ycol, header="", ybounds=None):
    """Emit speed_pct,<xcol>,<ycol> sorted by speed then x; skip non-numeric datasets.

    Points whose y falls outside ybounds are dropped as misclicks (e.g. the compressor 80%
    stray at PR=0.991) and reported.
    """
    rows, dropped = [], []
    for name, arr in datasets.items():
        spd = norm_speed(name)
        if spd is None:
            continue
        for x, y in arr[np.argsort(arr[:, 0])]:
            if ybounds is not None and not (ybounds[0] <= y <= ybounds[1]):
                dropped.append((spd, x, y))
                continue
            rows.append((spd, x, y))
    rows.sort(key=lambda r: (r[0], r[1]))
    with open(path, "w", newline="") as f:
        f.write(header)
        w = csv.writer(f)
        w.writerow(["speed_pct", xcol, ycol])
        for spd, x, y in rows:
            w.writerow([f"{spd:g}", f"{x:.5f}", f"{y:.5f}"])
    note = f"  (dropped {len(dropped)}: {[tuple(np.round(d,3)) for d in dropped]})" if dropped else ""
    print(f"  wrote {os.path.relpath(path, os.path.join(HERE, '..'))}  ({len(rows)} rows){note}")


def main():
    comp8 = parse_wpd(os.path.join(COMP_DIR, "fig 8.csv"))
    comp10 = parse_wpd(os.path.join(COMP_DIR, "fig 10.csv"))
    surge = parse_wpd(os.path.join(COMP_DIR, "fig 8 surge line.csv"))
    turb8 = parse_wpd(os.path.join(TURB_DIR, "fig 8.csv"))
    turb11a = parse_wpd(os.path.join(TURB_DIR, "fig 11a.csv"))
    turb11b = parse_wpd(os.path.join(TURB_DIR, "fig 11b.csv"))

    summarize("COMPRESSOR Fig 8 (PR vs flow)", comp8,
              "equiv flow W√θ/δ, lb/s", "PR P6/P1",
              design_spec=(0.581, 1.90), ybounds=(1.0, 2.2))
    summarize("COMPRESSOR Fig 10 (eff vs flow)", comp10,
              "equiv flow W√θ/δ, lb/s", "adiabatic eff",
              design_spec=(0.581, 0.795), ybounds=(0.0, 1.0))
    summarize("COMPRESSOR Fig 8 surge line", surge, "flow", "PR", ybounds=(1.0, 2.2))
    summarize("TURBINE Fig 8 (flow vs PR)", turb8,
              "PR (total-to-static)", "equiv flow, lb/s",
              design_spec=(1.69, 0.486), ybounds=(0.0, 1.0))
    summarize("TURBINE Fig 11a (static eff vs nu)", turb11a,
              "blade-jet speed ratio nu", "eta (total-to-static)",
              design_spec=(0.70, 0.872), ybounds=(0.0, 1.0))
    summarize("TURBINE Fig 11b (total eff vs nu)", turb11b,
              "blade-jet speed ratio nu", "eta (total-to-total)",
              design_spec=(0.70, 0.913), ybounds=(0.0, 1.0))

    print("\n---- writing authoritative digitized CSVs (data/digitized/) ----")
    h_cpr = ("# Compressor PR map, argon rig. NASA TM X-2129 (1970) Figure 8 (PR=P6/P1 vs equiv flow).\n"
             + PROV + "# speed_pct, We_lbps (equiv flow W*sqrt(theta)/delta, lb/s), PR\n")
    write_speedlines(comp8, os.path.join(DIG, "comp_pr_argon.csv"), "We_lbps", "PR",
                     header=h_cpr, ybounds=(1.0, 2.2))
    h_cef = ("# Compressor adiabatic efficiency map, argon rig. NASA TM X-2129 (1970) Figure 10.\n"
             + PROV + "# speed_pct, We_lbps (equiv flow, lb/s), eff (eta_1-6)\n")
    write_speedlines(comp10, os.path.join(DIG, "comp_eff_argon.csv"), "We_lbps", "eff",
                     header=h_cef, ybounds=(0.0, 1.0))
    h_tmf = ("# Turbine equivalent mass flow map, argon rig. NASA TN D-5090 (1969) Figure 8.\n"
             + PROV + "# Wc_lbps = eps*W*sqrt(theta_cr)/delta; PR_ts = inlet-total/exit-static.\n"
             "# speed_pct, PR_ts, Wc_lbps\n")
    write_speedlines(turb8, os.path.join(DIG, "turbine_massflow_argon.csv"), "PR_ts", "Wc_lbps",
                     header=h_tmf, ybounds=(0.0, 1.0))
    h_tts = ("# Turbine STATIC efficiency vs blade-jet speed ratio nu, argon rig.\n"
             "# NASA TN D-5090 (1969) Figure 11a (total-to-static). All speed lines collapse to one nu curve.\n"
             + PROV + "# speed_pct, nu, eta_ts\n")
    write_speedlines(turb11a, os.path.join(DIG, "turbine_eff_ts_nu_argon.csv"), "nu", "eta_ts",
                     header=h_tts, ybounds=(0.0, 1.0))
    h_ttt = ("# Turbine TOTAL efficiency vs blade-jet speed ratio nu, argon rig.\n"
             "# NASA TN D-5090 (1969) Figure 11b (total-to-total). All speed lines collapse to one nu curve.\n"
             + PROV + "# speed_pct, nu, eta_tt\n")
    write_speedlines(turb11b, os.path.join(DIG, "turbine_eff_tt_nu_argon.csv"), "nu", "eta_tt",
                     header=h_ttt, ybounds=(0.0, 1.0))
    # surge: single dataset
    sarr = list(surge.values())[0]
    with open(os.path.join(DIG, "comp_surge_argon.csv"), "w", newline="") as f:
        f.write("# Compressor surge line, argon rig. NASA TM X-2129 (1970) Figure 8 dashed line.\n" + PROV
                + "# We_lbps (equiv flow, lb/s), PR at surge\n")
        w = csv.writer(f); w.writerow(["We_lbps", "PR"])
        for x, y in sarr[np.argsort(sarr[:, 0])]:
            w.writerow([f"{x:.5f}", f"{y:.5f}"])
    print(f"  wrote data/digitized/comp_surge_argon.csv  ({len(sarr)} rows)")

    # ---- validation figure ----
    fig, ax = plt.subplots(2, 3, figsize=(18, 10))
    cmap = plt.cm.viridis

    def plot_lines(a, dsets, title, xl, yl, dspec=None):
        spds = sorted([norm_speed(n) for n in dsets if norm_speed(n) is not None])
        for name, arr in dsets.items():
            s = norm_speed(name)
            arr = arr[np.argsort(arr[:, 0])]
            if s is None:
                a.plot(arr[:, 0], arr[:, 1], "r*", ms=16, label=name, zorder=5)
            else:
                c = cmap((spds.index(s)) / max(1, len(spds) - 1))
                a.plot(arr[:, 0], arr[:, 1], "o-", color=c, ms=4, lw=1, label=f"{s:g}%")
        if dspec:
            a.plot(*dspec, "kx", ms=12, mew=2, label="spec")
        a.set(title=title, xlabel=xl, ylabel=yl); a.grid(alpha=0.3); a.legend(fontsize=7, ncol=2)

    plot_lines(ax[0, 0], comp8, "Comp Fig 8: PR vs flow", "flow lb/s", "PR", (0.581, 1.90))
    plot_lines(ax[0, 1], comp10, "Comp Fig 10: eff vs flow", "flow lb/s", "eff", (0.581, 0.795))
    # surge overlaid on comp8 panel
    ax[0, 0].plot(sarr[:, 0], sarr[:, 1], "k--", lw=1, label="surge")
    ax[0, 2].plot(sarr[:, 0], sarr[:, 1], "ko-", ms=3, lw=1); ax[0, 2].set(
        title="Comp Fig 8 surge line", xlabel="flow lb/s", ylabel="PR"); ax[0, 2].grid(alpha=0.3)
    plot_lines(ax[1, 0], turb8, "Turb Fig 8: flow vs PR", "PR (t-s)", "flow lb/s", (1.69, 0.486))
    plot_lines(ax[1, 1], turb11a, "Turb Fig 11a: eta_ts vs nu", "nu", "eta_ts", (0.70, 0.872))
    plot_lines(ax[1, 2], turb11b, "Turb Fig 11b: eta_tt vs nu", "nu", "eta_tt", (0.70, 0.913))
    fig.tight_layout()
    outpng = os.path.join(MD, "manual_digitize_check.png")
    fig.savefig(outpng, dpi=100); plt.close(fig)
    print(f"\nwrote {os.path.relpath(outpng, MD)}")


if __name__ == "__main__":
    main()
