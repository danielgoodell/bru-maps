"""
Reynolds-number and axial-clearance corrections for the BRU 4.97-in radial-inflow turbine.

REYNOLDS — NASA TN D-5090 (1969), Fig 14 (a faired cross-plot; "does not show any actual data
  points"). Definition (report):  Re = W / (mu * r_t)  [W mass flow, mu viscosity, r_t tip radius].
  Design Re = 76 200.  Total efficiency eta_t(1->3): 0.890 @ Re 34 950 -> 0.913 @ 76 200 ->
  0.925 @ 175 800 (a ~0.035 rise).  Loss power law (same form as the compressor):
      (1 - eta)/(1 - eta_ref) = (Re_ref/Re)^n,  Re_ref = 76 200, eta_ref = 0.913.
  PIECEWISE n (a single n cannot match both ends): the slope of ln(1-eta) vs ln(Re) steepens at
  low Re. Fitting the two segments through the three anchors (breakpoint at design Re):
      n = 0.30  for Re < 76 200   (low-Re, more loss-sensitive)
      n = 0.18  for Re >= 76 200  (high-Re, loss flattening toward fully-turbulent)
  These reproduce 0.890/0.913/0.925 exactly; the old single n = 0.22 gave 0.897/0.913/0.928.
  NOTE: static eff (1->3) is LESS Reynolds-sensitive (0.850->0.885, design 0.872 => n ~ 0.20 low /
  0.13 high). The map's effMap is TOTAL eff, so the total exponents below drive the map; the static
  exponents are recorded (N_RE_STATIC_*) for the tt-PR reconciliation if ever Reynolds-corrected.

CLEARANCE — NASA TM X-52552 (1969), the 5-in (12.62-cm) turbine = OURS.
  Total efficiency decrease ~ 1/3 percent (RELATIVE to design eta) per 1 percent increase in
  axial clearance (% of rotor blade height); mass-flow loss ~ 0.2 percent per 1 percent.
      eta  *= 1 - 0.00333*(clearance% - clearance_ref%)
      flow *= 1 - 0.00200*(clearance% - clearance_ref%)
  Cold-rig reference clearance ~ 2% (clearance report low point 1.94%); in-loop BRU ~ 12%
  (TM X-2350). Check: 2->12% => -3.33% relative => 0.894 -> 0.864 (TM X-2350 measured 0.860). OK.
"""
import numpy as np

R_T_M = 0.0631          # tip radius, m (4.97 in dia / 2)
RE_DES_TURB = 76200.0
ETA_REF_TURB = 0.913    # total eff (1->3) at design Re
N_RE_LOW = 0.30         # total-eff loss-law exponent for Re <  RE_DES_TURB
N_RE_HIGH = 0.18        # total-eff loss-law exponent for Re >= RE_DES_TURB
N_RE_STATIC_LOW = 0.20  # static-eff exponents (documented; map uses total above)
N_RE_STATIC_HIGH = 0.13

CLEAR_DETA_PER_PCT = 0.00333   # relative eff loss per 1% clearance increase (5-in turbine)
CLEAR_DFLOW_PER_PCT = 0.00200  # relative flow loss per 1% clearance increase
CLEAR_REF_PCT = 2.0            # cold-rig reference clearance (% blade height)


def reynolds_turbine(W_kgps, mu_pas, r_t=R_T_M):
    """Turbine Reynolds number Re = W/(mu*r_t)."""
    return W_kgps / (mu_pas * r_t)


def reynolds_loss_scale(Re, n_low=N_RE_LOW, n_high=N_RE_HIGH):
    """(1-eta)/(1-eta_ref) at Re relative to design Re (piecewise-n loss law, continuous at
    Re=RE_DES_TURB where the scale is 1.0). Pass the static exponents for static efficiency."""
    Re = np.asarray(Re, float)
    n = np.where(Re < RE_DES_TURB, n_low, n_high)
    return (RE_DES_TURB / Re) ** n


def correct_efficiency_reynolds(eta_map, Re):
    """Map efficiency (at design Re) -> efficiency at Re."""
    return 1.0 - (1.0 - np.asarray(eta_map)) * reynolds_loss_scale(Re)


def correct_clearance(eta_map, Wc_map, clearance_pct):
    """Apply axial-clearance penalty (eff and flow) relative to the cold-rig reference."""
    d = clearance_pct - CLEAR_REF_PCT
    eta = np.asarray(eta_map) * (1.0 - CLEAR_DETA_PER_PCT * d)
    wc = np.asarray(Wc_map) * (1.0 - CLEAR_DFLOW_PER_PCT * d)
    return eta, wc


if __name__ == "__main__":
    # Reynolds self-check (piecewise n now reproduces all three anchors exactly)
    for Re in (34950, 50000, 76200, 120000, 175800):
        eta = 1 - (1 - ETA_REF_TURB) * float(reynolds_loss_scale(Re))
        print(f"Re={Re:7d}  eta_t = {eta:.4f}")
    print("  (report anchors: 0.890 @ 34950, 0.913 @ 76200, 0.925 @ 175800)")
    # Clearance self-check vs TM X-2350 (apply to overall 1->4 design eff 0.894)
    eta12, _ = correct_clearance(0.894, 1.0, 12.0)
    print(f"\nclearance 2->12%: overall eta 0.894 -> {float(eta12):.3f}  (TM X-2350 measured 0.860)")
    eta12_13, w12 = correct_clearance(0.913, 1.0, 12.0)
    print(f"clearance 2->12%: map eta(1->3) 0.913 -> {float(eta12_13):.3f}, flow x{float(w12):.3f}")
