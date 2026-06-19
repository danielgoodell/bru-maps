"""
Emit NPSS-native Reynolds-number-index (RNI) correction tables for the BRU compressor and turbine.

NPSS applies Reynolds effects through the element's S_Re socket as DIRECT multipliers on the map
efficiency and corrected flow (verified against the NPSS docs and the T-MATS JT9D example maps):

    effBase = s_effDes * s_effRe * effMap          (compressor and turbine)
    WcBase  = s_WcDes  * s_WcRe  * WcMap           (compressor;  turbine uses WpMap / s_WpRe)

Only that usage is standard. NPSS specifies no convention for HOW the S_Re subelement produces the two
multipliers -- the internals are implementation-defined -- so we choose one: lookup tables returning the
multipliers from the Reynolds Number Index

    RNI = Re / ReDes        (ReDes is declared in the map; RNI = 1 at the reference Reynolds)

(normalized like the rest of the map, matches the one developer-guide example, and puts the design point
at RNI = 1 -> unit multipliers). Indexing on absolute Re instead is a trivial swap (the CSV carries both).
If the S_Re socket is left empty NPSS sets s_effRe = s_WcRe = 1.0 (no correction). NPSS corrects
EFFICIENCY and FLOW only -- there is no PR Reynolds scalar -- so the small compressor PR knock-down
the report notes is not represented here (it remains in the pre-baked `variants/loop_Re*` maps).

Both BRU reports give the effect in efficiency-DEFICIT form (1-eta)/(1-eta_ref) = (Re_ref/Re)^n =: S.
NPSS's s_effRe is a single multiplier on the whole effMap, so we anchor it at the loss law's own
reference efficiency eta_ref:

    s_effRe(RNI) = [1 - (1 - eta_ref) * S(RNI)] / eta_ref     (= 1.0 at RNI = 1)

exact at eta_ref and a close approximation elsewhere (NPSS's uniform-multiplier model). The deficit
multiplier S(RNI) is also written to the CSV so an S_Re subelement that reads the live effMap can be
exact instead:  s_effRe = [1 - (1 - effMap) * S(RNI)] / effMap.

Reference Reynolds numbers (ReDes, RNI = 1):
    compressor  Re_U = 3.03e6  (rho1*U_t3*D_t3/mu1; design ~3.1e6)  -- NASA TN D-6640
    turbine     Re   = 76 200  (W/(mu*r_t); total eff 1->3)          -- NASA TN D-5090

The S_Re subelement block is embedded directly in the base maps (the stage-1 builders import
compressor_s_re_block / turbine_s_re_block), so each map file is self-contained. This script, run on
its own, writes only the human-readable table dumps:
    maps/compressor/compressor_reynolds_correction.csv
    maps/turbine/turbine_reynolds_correction.csv
"""
import os
import numpy as np
import pandas as pd

import reynolds_correction as rc        # compressor (TN D-6640)
import turbine_corrections as tc        # turbine    (TN D-5090)

HERE = os.path.dirname(__file__)
CMAP = os.path.join(HERE, "..", "..", "maps", "compressor")
TMAP = os.path.join(HERE, "..", "..", "maps", "turbine")

ETA_REF_COMP = 1.0 - rc.LOSS_REF        # 0.813, the loss law's reference eta_max
ETA_REF_TURB = tc.ETA_REF_TURB          # 0.913


def _rni_grid(lo, hi, n, anchors):
    """Log-spaced RNI grid with the piecewise-law breakpoints forced in (kills interp error at kinks)."""
    g = np.geomspace(lo, hi, n)
    g = np.unique(np.round(np.concatenate([g, anchors]), 4))
    return g[(g >= lo) & (g <= hi)]


def _table1d(name, indep, x, y):
    """One NPSS 1-D table inside the S_Re subelement: name = f(indep)."""
    xs = ", ".join(f"{v:.4f}" for v in x)
    ys = ", ".join(f"{v:.5f}" for v in y)
    return (f"   Table {name}(real {indep}) {{\n"
            f"      {indep} = {{ {xs} }}\n"
            f"      {name} = {{ {ys} }}\n"
            f"   }}")


def _s_re_subelement(cls, tables):
    """Wrap the correction tables in the nested 'Subelement <cls> S_Re { ... }' block. Per the NPSS
    user manual the Reynolds subelement nests inside the map subelement (Subelement <MapType> S_map {
    Subelement <cls> S_Re { ... } }), filling its S_Re socket. The base-map builders embed this block
    directly so each map file is self-contained (build_compressor_map.py / build_turbine_map.py)."""
    return "\n".join([f"Subelement {cls} S_Re {{", *tables, "}"])


def _indent(block, n=3):
    """Indent every non-blank line of a text block by n spaces (to nest it inside S_map)."""
    pad = " " * n
    return "\n".join((pad + ln) if ln else ln for ln in block.split("\n"))


def compressor_s_re_block(indent=3):
    """The compressor Reynolds S_Re subelement, indented for embedding inside a CompressorRlineMap
    S_map block. Imported by build_compressor_map.py."""
    return _indent(compressor_table()[1], indent)


def turbine_s_re_block(indent=3):
    """The turbine Reynolds S_Re subelement, indented for embedding inside a TurbinePRmap S_map
    block (same block for the ts and tt maps -- Reynolds acts on efficiency only). Imported by
    build_turbine_map.py."""
    return _indent(turbine_table()[1], indent)


def _s_effRe(S, eta_ref):
    """Direct efficiency multiplier from the deficit scaling S, anchored at eta_ref. =1 at S=1."""
    return (1.0 - (1.0 - eta_ref) * S) / eta_ref


# ----------------------------------------------------------------------------- compressor
def compressor_table():
    edges = np.array([0.34e6, 0.43e6, 1.21e6, rc.RE_REF]) / rc.RE_REF
    RNI = _rni_grid(0.11, 2.5, 14, edges)
    Re = RNI * rc.RE_REF
    S = rc.loss_scale(Re)                       # efficiency-deficit scaling (1-eta)/(1-eta_ref)
    s_effRe = _s_effRe(S, ETA_REF_COMP)
    s_WcRe = np.ones_like(Re)                    # no Reynolds flow shift (surge/flow shift not applied)

    hdr = (
        "// NASA BRU 4.25-in centrifugal compressor -- Reynolds (RNI) correction for the S_Re socket.\n"
        "//   Socket: S_Re   Socket type: COMPRESSOR_REYNOLDS_EFFECTS   Returns: s_WcRe, s_effRe\n"
        "// Native form (per the NPSS user manual) is a subelement NESTED inside the map subelement:\n"
        "//   Subelement CompressorRlineMap S_map { ... Subelement CompressorReynoldsEffects S_Re {..} }\n"
        "//   -- paste the S_Re block below inside this map's S_map { } block. (Class name inferred from\n"
        "//   the socket type; the manual's worked example is the turbine's TurbineReynoldsEffects.)\n"
        "// NPSS fixes ONLY how the outputs are used: effBase = s_effDes*s_effRe*effMap, WcBase =\n"
        "//   s_WcDes*s_WcRe*WcMap. The subelement internals are implementation-defined (no NPSS\n"
        "//   standard). OUR convention: lookup tables returning the multipliers vs RNI = Re_U/ReDes,\n"
        "//   so RNI = 1 (s_effRe = s_WcRe = 1) at the reference Re. Re_U = rho1*U_t3*D_t3/mu1.\n"
        "// Source: NASA TN D-6640 (1972). NPSS applies these as DIRECT multipliers on the base map:\n"
        "//   effBase = s_effDes*s_effRe*effMap ;  WcBase = s_WcDes*s_WcRe*WcMap   (eff & flow only)\n"
        f"//   RNI = Re_U / ReDes,  ReDes = {rc.RE_REF:.3g}  (Re_U = rho1*U_t3*D_t3/mu1; design ~3.1e6).\n"
        f"// s_effRe is anchored at eta_ref = {ETA_REF_COMP:.3f}: = [1-(1-eta_ref)*S]/eta_ref, S the\n"
        "//   report's deficit law (1-eta)/(1-eta_ref)=(ReDes/Re)^n, piecewise n=0.06/0.09/0.20 as Re\n"
        f"//   falls. At the lowest tested Re (RNI~{0.34e6/rc.RE_REF:.2f}) s_effRe~{_s_effRe(rc.loss_scale(0.34e6).item(),ETA_REF_COMP):.3f}"
        f" -> eta_max {ETA_REF_COMP:.3f}*s_effRe = {ETA_REF_COMP*_s_effRe(rc.loss_scale(0.34e6).item(),ETA_REF_COMP):.3f}.\n"
        "// s_WcRe = 1 (report's surge/max-flow shift left for loop calibration). NPSS has no PR\n"
        "//   Reynolds scalar; the small peak-PR knock-down stays in the pre-baked variants/ maps.\n"
        "// (CSV carries the deficit S for an S_Re subelement that wants the exact per-point form.)\n"
    )
    body = _s_re_subelement("CompressorReynoldsEffects", [
        _table1d("s_effRe", "RNI", RNI, s_effRe),
        "",
        _table1d("s_WcRe", "RNI", RNI, s_WcRe),
    ])
    df = pd.DataFrame({"RNI": RNI, "Re_U": Re, "s_effRe": s_effRe, "s_WcRe": s_WcRe,
                       "lossDeficitS": S, "eta_max": ETA_REF_COMP * s_effRe})
    return hdr, body, df


# ----------------------------------------------------------------------------- turbine
def turbine_table():
    RNI = _rni_grid(0.30, 3.0, 12, np.array([1.0]))   # breakpoint at design Re (RNI=1)
    Re = RNI * tc.RE_DES_TURB
    S = tc.reynolds_loss_scale(Re)                    # total eff (1->3) deficit scaling
    s_effRe = _s_effRe(S, ETA_REF_TURB)
    s_WpRe = np.ones_like(Re)                          # Reynolds moves eff only, not flow (clearance does)

    hdr = (
        "// NASA BRU 4.97-in radial-inflow turbine -- Reynolds (RNI) correction for the S_Re socket.\n"
        "//   Socket: S_Re   Returns: s_WpRe, s_effRe  (turbine analog of COMPRESSOR_REYNOLDS_EFFECTS).\n"
        "// Native form (per the NPSS user manual) is a subelement NESTED inside the map subelement:\n"
        "//   Subelement TurbinePRmap S_map { ... Subelement TurbineReynoldsEffects S_Re { ... } }\n"
        "//   -- paste the S_Re block below inside this map's S_map { } block.\n"
        "// NPSS fixes ONLY how the outputs are used: effBase = s_effDes*s_effRe*effMap (and the WpMap\n"
        "//   analog). The subelement internals are implementation-defined (no NPSS standard). OUR\n"
        "//   convention: lookup tables returning the multipliers vs RNI = Re/ReDes, so RNI = 1\n"
        "//   (s_effRe = s_WpRe = 1) at design Re. Re = W/(mu*r_t).\n"
        "// Source: NASA TN D-5090 (1969) Fig 14. NPSS applies these as DIRECT multipliers:\n"
        "//   effBase = s_effDes*s_effRe*effMap ;  WpBase = s_WpDes*s_WpRe*WpMap\n"
        f"//   RNI = Re / ReDes,  ReDes = {tc.RE_DES_TURB:.0f}  (Re = W/(mu*r_t); RNI=1 at design Re).\n"
        f"// s_effRe anchored at eta_ref = {ETA_REF_TURB:.3f} (TOTAL eta 1->3): = [1-(1-eta_ref)*S]/eta_ref,\n"
        "//   S = (ReDes/Re)^n, piecewise n=0.30 (Re<ReDes)/0.18 (Re>=ReDes). Reproduces eta_t\n"
        "//   0.890/0.913/0.925 at Re 34950/76200/175800. (Static eff less Re-sensitive: n~0.20/0.13.)\n"
        "// s_WpRe = 1: Reynolds moves efficiency only; the flow penalty is axial CLEARANCE (separate).\n"
    )
    body = _s_re_subelement("TurbineReynoldsEffects", [
        _table1d("s_effRe", "RNI", RNI, s_effRe),
        "",
        _table1d("s_WpRe", "RNI", RNI, s_WpRe),
    ])
    df = pd.DataFrame({"RNI": RNI, "Re": Re, "s_effRe": s_effRe, "s_WpRe": s_WpRe,
                       "lossDeficitS": S, "eta_t": 1.0 - (1.0 - ETA_REF_TURB) * S})
    return hdr, body, df


def main():
    # The S_Re subelement is embedded directly in the base maps by build_{compressor,turbine}_map.py
    # (so each map file is self-contained). This script now only writes the human-readable CSVs of the
    # correction (s_effRe / s_WcRe / deficit vs RNI) for inspection and re-use by other tools.
    _, _, cdf = compressor_table()
    cdf.to_csv(os.path.join(CMAP, "compressor_reynolds_correction.csv"), index=False)
    _, _, tdf = turbine_table()
    tdf.to_csv(os.path.join(TMAP, "turbine_reynolds_correction.csv"), index=False)

    print("Wrote compressor_reynolds_correction.csv and turbine_reynolds_correction.csv")
    print("(the S_Re socket itself is embedded in the base maps by the map builders)")
    print(f"\nCompressor s_effRe vs RNI (eta_ref={ETA_REF_COMP:.3f}):")
    for rni in (0.112, 0.20, 0.40, 1.0, 1.5):
        S = float(rc.loss_scale(rni * rc.RE_REF))
        print(f"  RNI={rni:5.3f}  s_effRe={_s_effRe(S,ETA_REF_COMP):.4f}  eta_max={ETA_REF_COMP*_s_effRe(S,ETA_REF_COMP):.4f}")
    print(f"\nTurbine s_effRe vs RNI (eta_ref={ETA_REF_TURB:.3f}):")
    for rni in (0.459, 0.70, 1.0, 1.5, 2.307):
        S = float(tc.reynolds_loss_scale(rni * tc.RE_DES_TURB))
        print(f"  RNI={rni:5.3f}  s_effRe={_s_effRe(S,ETA_REF_TURB):.4f}  eta_t={1-(1-ETA_REF_TURB)*S:.4f}")


if __name__ == "__main__":
    main()
