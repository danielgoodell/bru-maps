"""
Build an NPSS-format compressor map from the digitized BRU compressor curves.

Pipeline:
  1. Load digitized PR(speed, We) and eff(speed, We) point sets (argon-referenced).
  2. Per speed line, parameterize from surge (beta=0) to max flow (beta=1) and resample
     PR, We, eff onto a common beta grid -> resolves the multivalued speed lines so each
     (Nc, beta) yields a unique (Wc, PR, eff).  [standard NPSS R-line/beta map structure]
  3. Normalize: NcMap = speed%/100,  WcMap = We / We_design.
  4. Emit:
       maps/compressor/compressor_argon.map           - NPSS Table syntax (PR, eff, Wc vs Nc, beta)
       maps/compressor/compressor_map_gridded.csv      - portable long-format grid (for pyCycle/t-MATS)
       maps/compressor/compressor_map_validation.png   - regenerated map overlaid on digitized points

Because all candidate fluids are gamma=5/3 monatomic, this normalized map is reusable for
He-Xe and Kr; only a Reynolds-number correction (Phase 3) differs between fluids.
"""
import os
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from corrected_params import COMP, COMP_SPEED_LINES

HERE = os.path.dirname(__file__)
DIG = os.path.join(HERE, "..", "..", "data", "digitized")
MAPS = os.path.join(HERE, "..", "..", "maps", "compressor")
os.makedirs(MAPS, exist_ok=True)

BETA = np.linspace(0.0, 1.0, 11)          # R-line grid, 0=surge ... 1=max flow
SPEEDS = [50, 60, 70, 80, 90, 100]        # percent design equivalent speed
We_DES = COMP.We_des_lbps                 # 0.581 lb/s


def load():
    pr = pd.read_csv(os.path.join(DIG, "comp_pr_argon.csv"), comment="#")
    eff = pd.read_csv(os.path.join(DIG, "comp_eff_argon.csv"), comment="#")
    return pr, eff


def align_eff_flow(We_e, We_p):
    """
    Map the efficiency curve's flow axis onto the PR curve's flow span for one speed line.

    PR (Fig 8) and efficiency (Fig 10) are measured at the SAME test points, so each speed line
    shares one surge and one max-flow. The two figures were digitized separately, so their flow
    extents drifted apart at the choke end (the surge ends already coincide). Rescaling the eff
    flow axis onto the (refined) PR span restores that shared range -- so the map carries one
    consistent (We, PR, eff) per beta and the eff curve no longer needs clamping/extrapolation.
    """
    return We_p[0] + (We_e - We_e[0]) / (We_e[-1] - We_e[0]) * (We_p[-1] - We_p[0])


def resample_speed_line(pr, eff, spd):
    """
    Return (We, PR, eff) resampled on the common BETA grid for one speed line.

    PR(We) and eff(We) are each fit with a shape-preserving PCHIP (piecewise-cubic Hermite):
    smooth (C1, no faceting), monotone-respecting, and -- unlike a quartic or smoothing spline --
    it does NOT overshoot above the data or invent spurious extrema (both of which distort the
    asymmetric efficiency peak and the flat PR tops). PR is read at refined accuracy (+/-0.01); the
    eff curve is a first-pass read, so PCHIP's faithfulness (no smoothing of real features) is
    appropriate -- the remaining limiter is the eff DATA, not the fit.
    """
    p = pr[pr.speed_pct == spd].sort_values("We_lbps")
    e = eff[eff.speed_pct == spd].sort_values("We_lbps")
    We_p, PR_p = p.We_lbps.to_numpy(), p.PR.to_numpy()
    We_e, EF_e = e.We_lbps.to_numpy(), e.eff.to_numpy()
    pr_fit = PchipInterpolator(We_p, PR_p)
    ef_fit = PchipInterpolator(align_eff_flow(We_e, We_p), EF_e)
    # beta is linear in flow from surge (0) to max flow (1) along the PR speed line
    We_g = We_p[0] + BETA * (We_p[-1] - We_p[0])
    return We_g, pr_fit(We_g), ef_fit(We_g)


def build_grid(pr, eff):
    Wc = np.zeros((len(SPEEDS), len(BETA)))
    PR = np.zeros_like(Wc)
    EF = np.zeros_like(Wc)
    for i, spd in enumerate(SPEEDS):
        We_g, PR_g, EFF_g = resample_speed_line(pr, eff, spd)
        Wc[i] = We_g / We_DES
        PR[i] = PR_g
        EF[i] = EFF_g
    Nc = np.array(SPEEDS) / 100.0
    return Nc, Wc, PR, EF


def write_npss_map(Nc, Wc, PR, EF, path):
    def table(name, Z, fmt="{:.4f}"):
        lines = [f"Table {name}(real NcMap, real betaMap) {{"]
        for i, nc in enumerate(Nc):
            lines.append(f"   NcMap = {nc:.3f} {{")
            betas = ", ".join(f"{b:.2f}" for b in BETA)
            vals = ", ".join(fmt.format(v) for v in Z[i])
            lines.append(f"      betaMap = {{ {betas} }}")
            lines.append(f"      {name.split('.')[-1]} = {{ {vals} }}")
            lines.append("   }")
        lines.append("}")
        return "\n".join(lines)

    # design-point location on the grid (Nc=1.0, We=We_DES)
    i100 = SPEEDS.index(100)
    beta_des = float(np.interp(1.0, Wc[i100], BETA))  # Wc=1.0 at design on 100% line

    header = f"""// ===================================================================
// NASA BRU 4.25-in sweptback-bladed centrifugal compressor map
// Generated from digitized NASA TM X-2129 (1970) argon rig data.
// Independents:  NcMap = N/sqrt(theta) normalized to 100% line
//                betaMap (R-line) 0 = surge ... 1 = max flow
// Outputs:       PRmap = P6/P1,  effMap = adiabatic eff,  WcMap = (W*sqrt(theta)/delta)/We_des
// Design point:  We_des = {We_DES} lb/s (argon-ref), PR = {COMP.PR_des}, eff = {COMP.eff_des}
//                lands at NcMap = 1.000, betaMap ~ {beta_des:.3f}
//
// NPSS scaler guidance (scaled-map CMP element): choose s_Nc, s_Wc, s_PR, s_eff so that
// at the cycle design point the element's corrected speed/flow/PR/eff equal the BRU design
// values above. For a different working fluid (He-Xe, Kr) the map is unchanged; apply the
// Reynolds-number multiplier from the 1972 Re report to effMap (and small WcMap) only.
// Surge line: betaMap = 0 column (see also data/digitized/comp_surge_argon.csv).
// Speed lines: shape-preserving PCHIP fits of PR(We) and eff(We) per line (smooth, no overshoot).
// PR refined to +/-0.01; eff is a first-pass read (re-digitize for production efficiency accuracy).
// ===================================================================
"""
    body = "\n\n".join([
        table("S_map.PRmap", PR),
        table("S_map.effMap", EF),
        table("S_map.WcMap", Wc),
    ])
    with open(path, "w") as f:
        f.write(header + "\n" + body + "\n")
    return beta_des


def write_csv(Nc, Wc, PR, EF, path):
    rows = []
    for i, nc in enumerate(Nc):
        for j, b in enumerate(BETA):
            rows.append((nc, b, Wc[i, j], PR[i, j], EF[i, j]))
    pd.DataFrame(rows, columns=["NcMap", "betaMap", "WcMap", "PRmap", "effMap"]).to_csv(
        path, index=False)


def validation_plot(pr, eff, Nc, Wc, PR, EF, path):
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.5))
    cmap = plt.cm.viridis(np.linspace(0, 1, len(SPEEDS)))
    for i, spd in enumerate(SPEEDS):
        c = cmap[i]
        We_grid = Wc[i] * We_DES
        ax[0].plot(We_grid, PR[i], "-", color=c, lw=1.5)
        d = pr[pr.speed_pct == spd]
        ax[0].plot(d.We_lbps, d.PR, "o", color=c, ms=4, label=f"{spd}%")
        ax[1].plot(We_grid, EF[i], "-", color=c, lw=1.5)
        de = eff[eff.speed_pct == spd]
        ax[1].plot(de.We_lbps, de.eff, "o", color=c, ms=4, label=f"{spd}%")
    # surge line — TM X-2129 Fig 8 draws it as a single straight dashed line (the surge
    # points are near-collinear), so fit one straight line rather than connecting points.
    s = pd.read_csv(os.path.join(DIG, "comp_surge_argon.csv"), comment="#")
    m, b = np.polyfit(s.We_lbps, s.PR, 1)
    we_s = np.array([s.We_lbps.min(), s.We_lbps.max()])
    ax[0].plot(we_s, m * we_s + b, "k--", lw=1, label="surge")
    ax[0].plot(We_DES, COMP.PR_des, "r*", ms=14, label="design")
    ax[0].set(xlabel="Equiv. weight flow W√θ/δ (argon), lb/s", ylabel="PR  P6/P1",
              title="Compressor PR map — regenerated (line) vs digitized (pts)")
    ax[1].plot(We_DES, COMP.eff_des, "r*", ms=14)
    ax[1].set(xlabel="Equiv. weight flow W√θ/δ (argon), lb/s", ylabel="adiabatic eff",
              title="Compressor efficiency map")
    for a in ax:
        a.grid(alpha=0.3); a.legend(fontsize=7, ncol=2)
    fig.tight_layout(); fig.savefig(path, dpi=110); plt.close(fig)


def main():
    pr, eff = load()
    Nc, Wc, PR, EF = build_grid(pr, eff)
    beta_des = write_npss_map(Nc, Wc, PR, EF, os.path.join(MAPS, "compressor_argon.map"))
    write_csv(Nc, Wc, PR, EF, os.path.join(MAPS, "compressor_map_gridded.csv"))
    validation_plot(pr, eff, Nc, Wc, PR, EF,
                    os.path.join(MAPS, "compressor_map_validation.png"))
    print("Wrote maps/compressor/compressor_argon.map, compressor_map_gridded.csv, validation.png")
    print(f"Design point lands at NcMap=1.000, betaMap={beta_des:.3f}")
    print(f"Grid: {len(SPEEDS)} speed lines x {len(BETA)} beta points")


if __name__ == "__main__":
    main()
