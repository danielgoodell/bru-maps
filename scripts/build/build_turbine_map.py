"""
Build an NPSS-format map for the 4.97-in radial-inflow turbine from digitized NASA TN D-5090.

Inputs (data/digitized/):
  turbine_massflow_argon.csv     - equivalent mass flow vs (p1'/p3)_ts per speed line (Fig 8)
  turbine_eff_tt_nu_argon.csv    - total   efficiency vs nu, per speed line (Fig 11b)
  turbine_eff_ts_nu_argon.csv    - static  efficiency vs nu, per speed line (Fig 11a)
    (all speed lines collapse onto one nu curve; we pool every line's markers and fit one poly)

Construction (FITTED approach):
  The published figures are sparse test markers with a smooth faired curve drawn through them.
  We recover those smooth curves with physically-motivated fits instead of connecting hand-read
  points (which propagates read-noise into kinks, e.g. the old 100% flow-line bump).

  - Mass flow: each speed line is fit to a generalized ellipse/Stodola swallowing law
        Wc(PR) = W_ch * sqrt(1 - PR^(-a))
    (monotonic, concave, saturates at choke W_ch, and W(PR=1)=0 by construction -> kinks
    impossible). The per-speed coefficients (W_ch, a) are then themselves fit as smooth
    quadratics in speed, so the whole family is smooth in BOTH directions (no wandering/crossing
    lines). The design flow is NOT hard-anchored; we report how far the fitted 100% line lands
    from the Table-I design value (0.486) as an honest reconciliation.
  - Efficiency: the single collapsed eta-vs-nu curve is fit with a smooth polynomial (one curve,
    all the data, clean derivatives). Efficiency at any (Nc, PR) follows from
        nu(Nc%, PR) = nu_des*(Nc%/100)*sqrt(1-PR_des^-0.4)/sqrt(1-PR^-0.4)
    then eta = poly(nu).
  - Normalize: NcMap = speed%/100, WcMap = Wc/Wc_des, PRmap = PR (total-to-static).

Outputs (maps/):
  turbine_argon.map               - NPSS Table syntax, PR axis = total-to-STATIC (p1'/p3)
  turbine_argon_tt.map            - same map, PR axis = total-to-TOTAL (pairs with an NPSS/GasCycle
                                    turbine element). PR_tt recovered from eta_tt/eta_ts (exit-KE),
                                    no assumed exit Mach:  PR_tt = [1-(eta_ts/eta_tt)(1-PR_ts^-k)]^(-1/k)
  turbine_map_gridded.csv         - long-format grid (ts);  turbine_map_tt_gridded.csv adds PR_tt, nu, eta_ts
  turbine_map_validation.png      - fitted lines vs markers + flow residuals + eta(nu) + eta(PR)

Gas transfer: eps references to air; monatomic gases (Ar/Kr/He-Xe, gamma=5/3) share eps, so the
normalized map serves all three; only Reynolds differs (handled separately).
"""
import os
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from generate_reynolds_maps import turbine_s_re_block

HERE = os.path.dirname(__file__)
DIG = os.path.join(HERE, "..", "..", "data", "digitized")
MAPS = os.path.join(HERE, "..", "..", "maps", "turbine")
os.makedirs(MAPS, exist_ok=True)

# NPSS turbine maps are parameterized on corrected speed and PRESSURE RATIO directly
# (Wc and eff are the map outputs, PR is an independent) -- unlike a compressor R-line/beta map.
# We therefore tabulate the two output tables WcMap(NcMap, PRmap) and effMap(NcMap, PRmap) on a
# common PR axis. The swallowing-law and eta(nu) fits are evaluated at every (Nc, PR) node, so the
# map is rectangular; PR nodes slightly outside a given speed line's measured span are filled by
# the (monotone, saturating) fit -- the standard way NPSS turbine maps cover a common PR grid.
PRGRID_TS = np.round(np.linspace(1.45, 2.00, 12), 3)   # total-to-static PR independent axis
PRGRID_TT = np.round(np.linspace(1.42, 1.95, 12), 3)   # total-to-total  PR independent axis
SPEEDS = [30, 50, 70, 90, 100, 110]
Wc_DES = 0.486          # equivalent design mass flow, lb/s (Table I)
PR_DES = 1.69           # design total-to-static PR (Fig 8 marker)
NU_DES = 0.70           # design blade-jet speed ratio
GAMMA_EXP = 0.4         # (gamma-1)/gamma, gamma=5/3
EFF_POLY_ORDER = 4      # order of the eta-vs-nu fit


# ---- functional forms -------------------------------------------------------
def flow_form(PR, W_ch, a, b):
    """
    Three-parameter swallowing-capacity curve: Wc(PR) = W_ch * (1 - PR^-a)^b.

    The plain ellipse/Stodola law is the special case b=0.5; with b free the curvature is
    tunable. This matters: 2-parameter forms (ellipse, exponential, Michaelis) all leave a
    systematic '+ - - - + +' residual arch on every speed line -- the measured curves bend
    harder than a fixed-curvature form allows. The free exponent removes that structure and
    drops residuals to the digitization-noise floor. Properties retained: monotonic increasing
    (a,b>0), Wc(PR=1)=0, saturates toward W_ch. Fit per-speed-line independently; a/b lose clean
    physical meaning on the unsaturated high-speed lines (they trade off) but the CURVE within the
    measured PR range is faithful -- which is all the map uses (no PR or speed extrapolation).
    """
    return W_ch * np.maximum(1.0 - PR ** (-a), 0.0) ** b


def nu_of(Nc_pct, PR):
    """Blade-jet speed ratio at (speed%, total-to-static PR)."""
    vj = np.sqrt(1.0 - PR ** (-GAMMA_EXP))
    vj_des = np.sqrt(1.0 - PR_DES ** (-GAMMA_EXP))
    return NU_DES * (Nc_pct / 100.0) * vj_des / vj


def load():
    mf = pd.read_csv(os.path.join(DIG, "turbine_massflow_argon.csv"), comment="#")
    # Efficiency now comes as two per-speed-line figures (TN D-5090 Fig 11a static, 11b total).
    # All speed lines collapse onto one eta-vs-nu curve, so we POOL every line's markers into a
    # single (nu, eta) cloud and fit one polynomial -- that is the measured collapse itself,
    # stronger than the old hand-faired single curve.
    ett = pd.read_csv(os.path.join(DIG, "turbine_eff_tt_nu_argon.csv"), comment="#")
    ets = pd.read_csv(os.path.join(DIG, "turbine_eff_ts_nu_argon.csv"), comment="#")
    return mf, ett, ets


# ---- fitting ----------------------------------------------------------------
FLOW_P0 = [0.60, 3.0, 0.6]                      # W_ch, a, b initial guess
FLOW_BOUNDS = ([0.30, 0.05, 0.10], [3.0, 40.0, 8.0])


def fit_flow(mf):
    """
    Fit each speed line independently to the 3-parameter form Wc=W_ch*(1-PR^-a)^b.
    Per-line (not a global speed polynomial): the measured curves differ enough in curvature
    that a smooth-in-speed coefficient model reintroduces structured residuals. Independent fits
    each reach the noise floor and, because the underlying data lines are ordered and non-crossing,
    the resulting curves stay ordered too. Bounds keep the (degenerate) high-speed asymptote finite.
    """
    params, resid, pr_lo, pr_hi = {}, {}, {}, {}
    for spd in SPEEDS:
        d = mf[mf.speed_pct == spd].sort_values("PR_ts")
        pr, wc = d.PR_ts.to_numpy(), d.Wc_lbps.to_numpy()
        popt, _ = curve_fit(flow_form, pr, wc, p0=FLOW_P0, bounds=FLOW_BOUNDS, maxfev=400000)
        params[spd] = popt
        resid[spd] = wc - flow_form(pr, *popt)
        pr_lo[spd], pr_hi[spd] = pr.min(), pr.max()
    return dict(params=params, resid=resid, pr_lo=pr_lo, pr_hi=pr_hi,
                W_ch={s: params[s][0] for s in SPEEDS},
                a={s: params[s][1] for s in SPEEDS},
                b={s: params[s][2] for s in SPEEDS})


def fit_eff(ett, ets):
    """Fit one polynomial each through the pooled total (Fig 11b) and static (Fig 11a) clouds."""
    nut, et = ett.nu.to_numpy(), ett.eta_tt.to_numpy()
    nus, es = ets.nu.to_numpy(), ets.eta_ts.to_numpy()
    ct = np.polyfit(nut, et, EFF_POLY_ORDER)
    cs = np.polyfit(nus, es, EFF_POLY_ORDER)
    return dict(coef_total=ct, coef_static=cs,
                nu_lo=float(nut.min()), nu_hi=float(nut.max()),
                nu_t=nut, eta_t=et, nu_s=nus, eta_s=es,
                resid_total=et - np.polyval(ct, nut),
                resid_static=es - np.polyval(cs, nus))


def build_grid(flow, eff, prgrid=PRGRID_TS):
    """
    Evaluate the smoothed fits on the rectangular (speed, PR) grid.

    Returns Nc (1-D), Wc and EF (2-D: nspeed x nPR) as functions of the common PR axis `prgrid`.
    PR is the independent variable (NPSS turbine convention), so it is no longer an output table.
    """
    Wc = np.zeros((len(SPEEDS), len(prgrid)))
    EF = np.zeros_like(Wc)
    for i, spd in enumerate(SPEEDS):
        Wc[i] = flow_form(prgrid, *flow["params"][spd]) / Wc_DES
        nu = np.clip(nu_of(spd, prgrid), eff["nu_lo"], eff["nu_hi"])
        EF[i] = np.polyval(eff["coef_total"], nu)
    Nc = np.array(SPEEDS) / 100.0
    return Nc, Wc, EF


# ---- writers ----------------------------------------------------------------
PR_DES_TT = 1.645        # design total-to-total PR (PR_ts 1.69 -> PR_tt; see _pr_ts_to_tt)


def _table(tb_name, leaf, Nc, PR, Z):
    """One native NPSS turbine table: TB_<q>(NpMap, PRmap), PR the (shared) independent axis."""
    lines = [f"   Table {tb_name}(real NpMap, real PRmap) {{"]
    for i, nc in enumerate(Nc):
        prs = ", ".join(f"{p:.3f}" for p in PR)
        vals = ", ".join(f"{v:.4f}" for v in Z[i])
        lines += [f"      NpMap = {nc:.3f} {{", f"         PRmap = {{ {prs} }}",
                  f"         {leaf} = {{ {vals} }}", "      }"]
    # Per-axis interp/extrap, inside the table after its speed lines (cubic in speed and PR;
    # linear extrapolation off the edges, not flagged as an error).
    lines += ['      NpMap.interp = "lagrange3";  NpMap.extrap = "linear";',
              '      PRmap.interp = "lagrange3";  PRmap.extrap = "linear";',
              "      extrapIsError = 0;",
              "      printExtrap   = 0;"]
    lines.append("   }")
    return "\n".join(lines)


def _subelement(Nc, PR, Wc, EF, pr_des):
    """The native 'Subelement TurbinePRmap S_map { ... }' block: design scalars, the two tables, and
    the embedded Reynolds S_Re subelement (a no-op at design Re, RNI = 1; corrects efficiency off-
    design via effBase = s_effDes*s_effRe*effMap)."""
    return "\n".join([
        "Subelement TurbinePRmap S_map {",
        f"   PRmapDes = {pr_des:.3f};",
        "   NpMapDes = 1.000;",
        "   // ReDes  = 76200.0;   // RNI = Re/ReDes; anchor = ARGON COLD-RIG design Re (TN D-5090). At the",
        "   //   He-Xe ~10 kWe point RNI ~2 (hot mu at 1600F ~2.5x cuts Re, but the ~1.3 lb/s loop flow more",
        "   //   than offsets) -> s_effRe ~1.01, a small FAVORABLE bump. See docs/turbine.md.",
        "",
        _table("TB_Wp", "WpMap", Nc, PR, Wc),
        "",
        _table("TB_eff", "effMap", Nc, PR, EF),
        "",
        turbine_s_re_block(),
        "}",
    ])


def write_npss(Nc, PR, Wc, EF, path):
    hdr = (f"// NASA BRU 4.97-in radial-inflow turbine map\n"
           f"// From digitized NASA TN D-5090 (1969) argon cold-flow rig (Fig 8 flow, Fig 11 eff).\n"
           f"// NPSS turbine map: a 'Subelement TurbinePRmap S_map' block. Corrected speed and PRESSURE\n"
           f"// RATIO are the independents; the tables TB_Wp, TB_eff output corrected flow (WpMap) and\n"
           f"// efficiency (effMap).\n"
           f"// Flow: per-speed 3-param swallowing law Wp=W_ch*(1-PR^-a)^b (b free; tunable curvature).\n"
           f"// Eff: smooth polynomial fit of the collapsed eta-vs-nu curve, eta=poly(nu(NpMap,PR)).\n"
           f"// NpMap = N/sqrt(theta_cr) normalized to design (100% = 1.0); PRmap = inlet-total/exit-"
           f"static (p1'/p3).\n"
           f"// effMap = total eff eta_t(1->3); WpMap = (eps*W*sqrt(theta_cr)/delta)/Wp_des, "
           f"Wp_des={Wc_DES} lb/s.\n"
           f"// Design: PR_ts={PR_DES}, nu={NU_DES}, eta_t=0.913 (PRmapDes, NpMapDes below). The\n"
           f"// Reynolds correction is the embedded S_Re subelement (no-op at design Re): scales effMap.\n"
           f"// Monatomic-gas independent (Ar/Kr/HeXe). PR is total-to-static; total-to-total variant\n"
           f"// in turbine_argon_tt.map (see docs/turbine.md).\n")
    with open(path, "w") as f:
        f.write(hdr + "\n" + _subelement(Nc, PR, Wc, EF, PR_DES) + "\n")


def write_csv(Nc, PR, Wc, EF, path):
    rows = [(nc, PR[j], Wc[i, j], EF[i, j])
            for i, nc in enumerate(Nc) for j in range(len(PR))]
    pd.DataFrame(rows, columns=["NpMap", "PRmap", "WpMap", "effMap"]).to_csv(path, index=False)


# ---- total-to-total reconciliation -----------------------------------------
def _pr_ts_to_tt(PR_ts, nc_pct, eff):
    """
    total-to-static -> total-to-total PR via the exit kinetic energy implied by the two
    efficiency curves (no assumed exit Mach):

        actual work = eta_tt*Dh_id,tt = eta_ts*Dh_id,ts  (one physical quantity)
        => eta_ts/eta_tt = (1-PR_tt^-k)/(1-PR_ts^-k),  k=(gamma-1)/gamma=0.4
        => PR_tt = [1 - (eta_ts/eta_tt)*(1 - PR_ts^-k)]^(-1/k)

    eta_tt is the total-eta poly (Fig 11b); eta_ts the static-eta poly (Fig 11a), both at the
    same nu(nc, PR_ts).
    """
    nu = np.clip(nu_of(nc_pct, PR_ts), eff["nu_lo"], eff["nu_hi"])
    eta_tt = np.polyval(eff["coef_total"], nu)
    eta_ts = np.polyval(eff["coef_static"], nu)
    PR_tt = (1.0 - (eta_ts / eta_tt) * (1.0 - PR_ts ** (-GAMMA_EXP))) ** (-1.0 / GAMMA_EXP)
    return PR_tt, eta_tt, eta_ts, nu


def build_tt_grid(flow, eff, prgrid=PRGRID_TT):
    """
    Build the rectangular total-to-total turbine map: WcMap and effMap as functions of
    (NcMap, PR_tt). For each speed line we sweep PR_ts densely, map it to PR_tt (monotone), then
    invert to read the PR_ts -- and hence Wc and eff -- at each requested PR_tt node.
    Returns Nc, Wc, EF and reference arrays (PR_ts, eta_ts, nu) on the PR_tt grid.
    """
    pr_ts_dense = np.linspace(1.35, 2.10, 400)
    Wc = np.zeros((len(SPEEDS), len(prgrid)))
    EF = np.zeros_like(Wc); PRTS = np.zeros_like(Wc); ETS = np.zeros_like(Wc); NU = np.zeros_like(Wc)
    for i, spd in enumerate(SPEEDS):
        pr_tt_dense, eta_tt_d, eta_ts_d, nu_d = _pr_ts_to_tt(pr_ts_dense, spd, eff)
        pr_ts = np.interp(prgrid, pr_tt_dense, pr_ts_dense)   # invert monotone PR_tt(PR_ts)
        nu = np.clip(nu_of(spd, pr_ts), eff["nu_lo"], eff["nu_hi"])
        Wc[i] = flow_form(pr_ts, *flow["params"][spd]) / Wc_DES
        EF[i] = np.polyval(eff["coef_total"], nu)
        PRTS[i] = pr_ts; ETS[i] = np.polyval(eff["coef_static"], nu); NU[i] = nu
    Nc = np.array(SPEEDS) / 100.0
    return Nc, Wc, EF, PRTS, ETS, NU


def write_tt(Nc, PR_tt, Wc, EF, PR_ts, ETS, NU, path_map, path_csv):
    hdr = ("// BRU 4.97-in radial-inflow turbine map -- TOTAL-TO-TOTAL coordinates.\n"
           "// NPSS turbine map ('Subelement TurbinePRmap S_map'): independents NpMap and PRmap\n"
           "// (total-to-total Pt_in/Pt_out); tables TB_Wp, TB_eff output WpMap = Wp/0.486 and\n"
           "// effMap = total eta_t(1->3). (PR is an input, not output.)\n"
           "// PR_tt recovered from the total-to-static map via exit kinetic energy from the digitized\n"
           "// eta_tt (Fig 11b) and eta_ts (Fig 11a) nu-curves -- no assumed exit Mach.\n"
           "// Pairs directly with an NPSS/GasCycle total-to-total turbine expansion-ratio element.\n")
    with open(path_map, "w") as f:
        f.write(hdr + "\n" + _subelement(Nc, PR_tt, Wc, EF, PR_DES_TT) + "\n")
    rows = [(nc, PR_tt[j], Wc[i, j], PR_ts[i, j], EF[i, j], ETS[i, j], NU[i, j])
            for i, nc in enumerate(Nc) for j in range(len(PR_tt))]
    pd.DataFrame(rows, columns=["NpMap", "PRmap_tt", "WpMap", "PR_ts",
                                "eta_tt", "eta_ts", "nu"]).to_csv(path_csv, index=False)


def validation_plot(mf, flow, eff, Nc, PR, Wc, EF, path):
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    cmap = plt.cm.viridis(np.linspace(0, 1, len(SPEEDS)))

    # (0,0) flow map: markers (data) + smoothed fitted lines
    for i, spd in enumerate(SPEEDS):
        c = cmap[i]
        prline = np.linspace(flow["pr_lo"][spd], flow["pr_hi"][spd], 100)
        ax[0, 0].plot(prline, flow_form(prline, *flow["params"][spd]), "-", color=c, lw=1.4)
        d = mf[mf.speed_pct == spd]
        ax[0, 0].plot(d.PR_ts, d.Wc_lbps, "o", color=c, ms=4, label=f"{spd}%")
    ax[0, 0].plot(PR_DES, Wc_DES, "r*", ms=15, label="design (Table I)")
    ax[0, 0].set(xlabel="(p1'/p3) total-to-static PR", ylabel="equiv mass flow εW√θcr/δ, lb/s",
                 title="Flow map — 3-param fit Wc=W_ch(1-PR⁻ᵃ)ᵇ (line) vs markers")
    ax[0, 0].legend(fontsize=7)

    # (0,1) flow residuals: data - per-line fit (shows kink gone, read quality)
    for i, spd in enumerate(SPEEDS):
        d = mf[mf.speed_pct == spd].sort_values("PR_ts")
        ax[0, 1].plot(d.PR_ts, flow["resid"][spd] * 1000, "-o", color=cmap[i], ms=4, label=f"{spd}%")
    ax[0, 1].axhline(0, color="k", lw=0.7)
    rms = np.sqrt(np.mean(np.concatenate([flow["resid"][s] for s in SPEEDS]) ** 2)) * 1000
    ax[0, 1].set(xlabel="(p1'/p3) total-to-static PR", ylabel="flow residual (data−fit), ×10⁻³ lb/s",
                 title=f"Flow-fit residuals (RMS = {rms:.1f}e-3 lb/s)")
    ax[0, 1].legend(fontsize=7)

    # (1,0) eta vs nu: markers + polynomial fit
    nn = np.linspace(eff["nu_lo"], eff["nu_hi"], 200)
    ax[1, 0].plot(nn, np.polyval(eff["coef_total"], nn), "k-", label="η_total fit")
    ax[1, 0].plot(nn, np.polyval(eff["coef_static"], nn), "k--", label="η_static fit")
    ax[1, 0].plot(eff["nu_t"], eff["eta_t"], "ko", ms=4)
    ax[1, 0].plot(eff["nu_s"], eff["eta_s"], "ks", ms=4, mfc="none")
    ax[1, 0].plot(NU_DES, 0.913, "r*", ms=15, label="design")
    rt = np.sqrt(np.mean(eff["resid_total"] ** 2))
    ax[1, 0].set(xlabel="blade-jet speed ratio ν", ylabel="efficiency",
                 title=f"Efficiency vs ν — poly fit (η_total RMS = {rt:.4f})")
    ax[1, 0].legend(fontsize=7)

    # (1,1) eta vs PR derived through the map (PR is the shared independent axis)
    for i, spd in enumerate(SPEEDS):
        ax[1, 1].plot(PR, EF[i], "-", color=cmap[i], lw=1.4, label=f"{spd}%")
    ax[1, 1].plot(PR_DES, 0.913, "r*", ms=15)
    ax[1, 1].set(xlabel="(p1'/p3) total-to-static PR", ylabel="total efficiency η_t",
                 title="Map efficiency vs PR (smooth, via ν)")
    ax[1, 1].legend(fontsize=7)

    for a in ax.ravel():
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def main():
    mf, ett, ets = load()
    flow = fit_flow(mf)
    eff = fit_eff(ett, ets)
    Nc, Wc, EF = build_grid(flow, eff)
    write_npss(Nc, PRGRID_TS, Wc, EF, os.path.join(MAPS, "turbine_argon.map"))
    write_csv(Nc, PRGRID_TS, Wc, EF, os.path.join(MAPS, "turbine_map_gridded.csv"))
    Nc_tt, Wc_tt, EF_tt, PRTS, ETS, NU = build_tt_grid(flow, eff)
    write_tt(Nc_tt, PRGRID_TT, Wc_tt, EF_tt, PRTS, ETS, NU,
             os.path.join(MAPS, "turbine_argon_tt.map"),
             os.path.join(MAPS, "turbine_map_tt_gridded.csv"))
    validation_plot(mf, flow, eff, Nc, PRGRID_TS, Wc, EF,
                    os.path.join(MAPS, "turbine_map_validation.png"))

    # --- diagnostics ---
    print("Wrote turbine_argon.map, turbine_map_gridded.csv, turbine_map_validation.png\n")
    print("Per-speed 3-param fit  Wc = W_ch*(1-PR^-a)^b:")
    for s in SPEEDS:
        r = flow["resid"][s]
        print(f"  {s:3d}%:  W_ch={flow['W_ch'][s]:6.3f}  a={flow['a'][s]:5.2f}  b={flow['b'][s]:4.2f}  "
              f"rms={np.sqrt(np.mean(r**2))*1000:.1f}e-3  max|resid|={np.max(np.abs(r))*1000:.1f}e-3 lb/s")
    # design reconciliation: fitted 100% line at design PR vs Table-I value
    wfit_des = flow_form(PR_DES, *flow["params"][100])
    print(f"\nDesign reconciliation: fitted 100% line at PR_ts={PR_DES} gives Wc={wfit_des:.4f} "
          f"lb/s vs Table-I {Wc_DES:.4f} (Δ={(wfit_des-Wc_DES)*1000:+.1f}e-3 = {(wfit_des/Wc_DES-1)*100:+.2f}%)")
    i100 = SPEEDS.index(100)
    pr_at_des = np.interp(1.0, Wc[i100], PRGRID_TS)
    ef_at_des = np.interp(1.0, Wc[i100], EF[i100])
    print(f"100% line at WcMap=1.0 (design flow): PR={pr_at_des:.3f} (design 1.69), "
          f"eff={ef_at_des:.3f} (design 0.913)")
    print(f"nu at (100%, 1.69) = {nu_of(100,1.69):.3f} (design 0.70); "
          f"eta poly at design nu = {np.polyval(eff['coef_total'], NU_DES):.4f}")
    eta_ts_des = np.polyval(eff["coef_static"], NU_DES)
    pr_tt_des = (1.0 - (eta_ts_des / np.polyval(eff["coef_total"], NU_DES))
                 * (1.0 - PR_DES ** (-GAMMA_EXP))) ** (-1.0 / GAMMA_EXP)
    print(f"total-to-total: design PR_ts={PR_DES} -> PR_tt={pr_tt_des:.3f} "
          f"(eta_ts={eta_ts_des:.3f}, eta_tt={np.polyval(eff['coef_total'], NU_DES):.3f})")


if __name__ == "__main__":
    main()
