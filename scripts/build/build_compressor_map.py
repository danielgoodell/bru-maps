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
from generate_reynolds_maps import compressor_s_re_block

HERE = os.path.dirname(__file__)
DIG = os.path.join(HERE, "..", "..", "data", "digitized")
MAPS = os.path.join(HERE, "..", "..", "maps", "compressor")
os.makedirs(MAPS, exist_ok=True)

BETA = np.linspace(0.0, 1.0, 11)          # internal resample parameter, 0=surge ... 1=max flow
# NPSS map R-line axis. NPSS compressor maps (CompressorRlineMap subelement "S_map") use an
# RlineMap independent that conventionally runs ~1.0 (surge) to ~3.0+ (choke), with the design
# operating point near RlineMap=2.0. We map the internal surge->max-flow parameter linearly onto
# that range so the file is a drop-in NPSS R-line map (surge column = RlineMap 1.0).
RLINE_LO, RLINE_HI = 1.0, 3.0
RLINE = RLINE_LO + BETA * (RLINE_HI - RLINE_LO)
# Surge/stall line. In an NPSS R-line map the surge line is not a separate curve: it is the
# lowest R-line (lowest corrected flow at each speed). NPSS reads the scalar S_map.RlineStall as
# the stall boundary and computes the stall margins SMW/SMN from it; pyCycle's CompressorRlineMap
# uses the same convention (RlineStall = 1.0, min RlineMap = 1.0). Our surge column is RLINE_LO.
RLINE_STALL = RLINE_LO
# Variable-geometry (alpha) axis. The BRU is fixed-geometry, so the map is alpha-independent, but
# NPSS still interpolates the alpha axis. We therefore emit TWO identical sub-maps (alphaMap = 0 and
# 1) with LINEAR alpha interpolation: whatever alpha the element passes returns the same map, so the
# model never has to pin alphaMap to 0. (Set ALPHA_LEVELS = [0.0] for a single sub-map instead.)
ALPHA_LEVELS = [0.0, 1.0]
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
    # One native NPSS map table: TB_<q>(alphaMap, NcMap, RlineMap) with the leaf value array named
    # for the quantity (WcMap / PRmap / effMap). alphaMap is the variable-geometry axis; the BRU is
    # fixed-geometry, so a single alphaMap = 0.0 sub-map carries the whole map (the standard way to
    # express a fixed-geometry compressor in the CompressorRlineMap form).
    def table(tb_name, leaf, Z, fmt="{:.4f}"):
        lines = [f"   Table {tb_name}(real alphaMap, real NcMap, real RlineMap) {{"]
        for a in ALPHA_LEVELS:
            lines.append(f"      alphaMap = {a:.3f} {{")
            for i, nc in enumerate(Nc):
                rlines = ", ".join(f"{r:.2f}" for r in RLINE)
                vals = ", ".join(fmt.format(v) for v in Z[i])
                lines += [f"         NcMap = {nc:.3f} {{",
                          f"            RlineMap = {{ {rlines} }}",
                          f"            {leaf} = {{ {vals} }}",
                          "         }"]
            lines.append("      }")
        # Per-axis interp/extrap methods live INSIDE each table, after its alpha sub-maps (each
        # table sets them independently): cubic in speed and R-line, linear in alpha; linear
        # extrapolation off the edges, not flagged as an error.
        lines += [
            '      alphaMap.interp = "linear";     alphaMap.extrap = "linear";',
            '      NcMap.interp    = "lagrange3";  NcMap.extrap    = "linear";',
            '      RlineMap.interp = "lagrange3";  RlineMap.extrap = "linear";',
            "      extrapIsError = 0;",
            "      printExtrap   = 0;",
        ]
        lines.append("   }")
        return "\n".join(lines)

    # design-point location on the grid (Nc=1.0, We=We_DES)
    i100 = SPEEDS.index(100)
    beta_des = float(np.interp(1.0, Wc[i100], BETA))  # Wc=1.0 at design on 100% line
    rline_des = RLINE_LO + beta_des * (RLINE_HI - RLINE_LO)

    header = f"""// ===================================================================
// NASA BRU 4.25-in sweptback-bladed centrifugal compressor map
// Generated from digitized NASA TM X-2129 (1970) argon rig data.
// NPSS R-line map: a "Subelement CompressorRlineMap S_map" block holding the three native tables
//   TB_Wc, TB_PR, TB_eff -- corrected flow, pressure ratio and efficiency as functions of the
//   variable-geometry angle, corrected speed and R-line.
// Independents:  alphaMap = variable-geometry angle (fixed-geometry BRU: a single 0.0 sub-map)
//                NcMap    = N/sqrt(theta) normalized to the 100% line
//                RlineMap = compressor R-line, 1.0 (surge) ... 3.0 (max flow)
//   (fixed-geometry BRU: two identical alphaMap = 0/1 sub-maps + linear alpha interp, so any
//    alpha returns the same map; per-axis interp/extrap methods are declared after the tables.)
// Leaf outputs:  WcMap = (W*sqrt(theta)/delta)/We_des,  PRmap = P6/P1,  effMap = adiabatic eff
// Design point:  We_des = {We_DES} lb/s (argon-ref), PR = {COMP.PR_des}, eff = {COMP.eff_des}
//                lands at NcMap = 1.000, RlineMap ~ {rline_des:.3f}  (RlineMapDes, below)
//
// NPSS scaler guidance (scaled-map Compressor element): choose s_Nc, s_Wc, s_PR, s_eff so that
// at the cycle design point the element's corrected speed/flow/PR/eff equal the BRU design
// values above, with RlineMapDes set to the value below. For a different working fluid (He-Xe,
// Kr) the map is unchanged; the Reynolds correction is the embedded S_Re subelement below
// (a no-op at the reference Re, RNI = 1) -- it scales effMap and WcMap only.
// Surge line: the RlineMap = {RLINE_STALL:.2f} column is the surge/stall boundary (each speed line's
// lowest corrected flow). RlineStall (below) tells NPSS where it is so the element computes the
// stall margins SMW (const flow) and SMN (const speed); the operating RlineMap is held >= it.
// The surge column traces the digitized Fig-8 surge line (data/digitized/comp_surge_argon.csv) to
// within read scatter (it is the last stable rig point, ~0.01-0.06 in PR inside the faired line).
// Speed lines: shape-preserving PCHIP fits of PR(We) and eff(We) per line (smooth, no overshoot).
// PR refined to +/-0.01; eff is a first-pass read (re-digitize for production efficiency accuracy).
// ===================================================================
"""
    # Native form: the design scalars are member declarations INSIDE the Subelement block (the way
    # the T-MATS example .map files carry them), not free S_map.<x> assignments.
    sub = [
        "Subelement CompressorRlineMap S_map {",
        "   alphaMapDes = 0.000;",
        "   NcMapDes    = 1.000;",
        f"   RlineMapDes = {rline_des:.3f};",
        f"   RlineStall  = {RLINE_STALL:.3f};",
        "   // ReDes    = 3030000.0;   // RNI = Re_U/ReDes; anchor = ARGON-RIG design Re_U (TN D-6640).",
        "   //   He-Xe operating points: RNI ~0.92 (6 kWe ref design) .. ~1.6 (10 kWe net, the real",
        "   //   operational design) -- the rig Re sits at the LOW end of that range. See docs/compressor.md.",
        "",
        table("TB_Wc", "WcMap", Wc),
        "",
        table("TB_PR", "PRmap", PR),
        "",
        table("TB_eff", "effMap", EF),
        "",
        "   // Reynolds correction filling the S_Re socket (RNI = Re_U/ReDes). A no-op at the argon-rig",
        "   // reference Re (RNI = 1); at the He-Xe ~10 kWe operational point RNI ~ 1.6 -> s_effRe ~ 1.006",
        "   // (a small FAVORABLE bump). NPSS applies effBase = s_effDes*s_effRe*effMap, WcBase similarly.",
        compressor_s_re_block(),
        "}",
    ]
    with open(path, "w") as f:
        f.write(header + "\n" + "\n".join(sub) + "\n")
    return rline_des


def write_csv(Nc, Wc, PR, EF, path):
    rows = []
    for i, nc in enumerate(Nc):
        for j, r in enumerate(RLINE):
            rows.append((nc, r, Wc[i, j], PR[i, j], EF[i, j]))
    pd.DataFrame(rows, columns=["NcMap", "RlineMap", "WcMap", "PRmap", "effMap"]).to_csv(
        path, index=False)


def verify_surge(Nc, Wc, PR):
    """
    Check the surge column (RlineMap = RLINE_STALL) against the independently digitized Fig-8
    surge line. Each speed line's beta=0 point should sit on that faired line; we report the PR
    residual (data minus faired-line PR at the same flow) so regenerations catch any drift.
    """
    s = pd.read_csv(os.path.join(DIG, "comp_surge_argon.csv"), comment="#").sort_values("We_lbps")
    print(f"\nSurge column (RlineMap={RLINE_STALL:.2f}) vs digitized Fig-8 surge line:")
    print(f"  {'Nc%':>4} {'We_surge':>9} {'PR_col':>8} {'PR_faired':>10} {'dPR':>8}")
    resid = []
    for i, nc in enumerate(Nc):
        we = Wc[i, 0] * We_DES                       # flow at the surge column
        pr_col = PR[i, 0]                            # speed line's surge PR
        pr_faired = float(np.interp(we, s.We_lbps, s.PR))
        resid.append(pr_col - pr_faired)
        print(f"  {nc*100:4.0f} {we:9.4f} {pr_col:8.4f} {pr_faired:10.4f} {pr_col-pr_faired:+8.4f}")
    print(f"  RMS dPR = {np.sqrt(np.mean(np.square(resid))):.4f}, "
          f"max |dPR| = {np.max(np.abs(resid)):.4f}")


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
    rline_des = write_npss_map(Nc, Wc, PR, EF, os.path.join(MAPS, "compressor_argon.map"))
    write_csv(Nc, Wc, PR, EF, os.path.join(MAPS, "compressor_map_gridded.csv"))
    validation_plot(pr, eff, Nc, Wc, PR, EF,
                    os.path.join(MAPS, "compressor_map_validation.png"))
    print("Wrote maps/compressor/compressor_argon.map, compressor_map_gridded.csv, validation.png")
    print(f"Design point lands at NcMap=1.000, RlineMap={rline_des:.3f} (set S_map.RlineMapDes)")
    print(f"Grid: {len(SPEEDS)} speed lines x {len(RLINE)} R-line points "
          f"({RLINE_LO:.1f}=surge/RlineStall .. {RLINE_HI:.1f}=max flow)")
    verify_surge(Nc, Wc, PR)


if __name__ == "__main__":
    main()
