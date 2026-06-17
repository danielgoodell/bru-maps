"""
Corrected / equivalent parameter conventions for the NASA BRU turbomachinery maps.

All NASA "equivalent" quantities in these reports use standard-day referral:
    T_std = 518.67 R, P_std = 14.696 psia
    theta = T1 / T_std,  delta = P1 / P_std
    W_corr = W * sqrt(theta) / delta      (equivalent weight flow)
    N_corr = N / sqrt(theta)              (equivalent speed)
    U_corr = U / sqrt(theta)              (equivalent tip speed)

Verified against TM X-2129 Table II (argon compressor design point):
    W=0.785 lb/s @ 20.25 psia / 540 R  ->  W_corr = 0.581 lb/s  (report: 0.581)

Gas note: argon, krypton, and the design He-Xe mixture are all monatomic (gamma = 5/3).
Argon was the rig surrogate chosen to match the design specific-heat ratio. Hence the
*normalized* corrected map is reusable across these fluids; only Reynolds number differs
(handled separately from the 1972 Reynolds-number report).
"""
from dataclasses import dataclass

T_STD_R = 518.67      # standard-day total temperature, deg R
P_STD_PSIA = 14.696   # standard-day total pressure, psia
GAMMA_MONATOMIC = 5.0 / 3.0

# Universal gas constant / molecular weights (lb-mol basis), for gas-property work
MW = {"argon": 39.948, "krypton": 83.80, "hexe_83p8": 83.8}


def theta(T1_R: float) -> float:
    return T1_R / T_STD_R


def delta(P1_psia: float) -> float:
    return P1_psia / P_STD_PSIA


def w_corr(W, T1_R, P1_psia):
    """Equivalent (corrected) weight flow W*sqrt(theta)/delta."""
    return W * theta(T1_R) ** 0.5 / delta(P1_psia)


def n_corr(N, T1_R):
    """Equivalent (corrected) speed N/sqrt(theta)."""
    return N / theta(T1_R) ** 0.5


@dataclass(frozen=True)
class CompressorDesign:
    """4.25-in sweptback-bladed centrifugal compressor design point (TM X-2129)."""
    We_des_lbps: float = 0.581      # equivalent weight flow at design (argon-referenced)
    N_des_rpm_equiv: float = 51176  # equivalent speed of the 100% line
    U_tip_des_equiv_fps: float = 949
    PR_des: float = 1.90
    eff_des: float = 0.795
    eff_peak_des_speed: float = 0.819   # at We = 0.508
    T6_T1_des: float = 1.37
    work_factor: float = 0.658
    Re_des: float = 3.8e6
    impeller_exit_dia_in: float = 4.25
    n_impeller_blades: int = 15
    backsweep_deg: float = 30.0
    n_diffuser_vanes: int = 17


COMP = CompressorDesign()

# Map speed-line key (TM X-2129 figures): U_tip,eq [ft/s] -> percent of design speed
COMP_SPEED_LINES = {100: 949, 90: 854, 80: 759, 70: 664, 60: 569, 50: 475}

# Gas-transfer scalars. The normalized compressor map (NcMap=N/N_des, WcMap=We/We_des) is
# gas-INDEPENDENT; a fluid enters only through its design equivalent flow & speed, which scale
# as sqrt(MW) relative to argon (verified vs the He-Xe loop map, NASA TM X-67989). Use these to
# set the NPSS map scalers (s_Wc, s_Nc) for the chosen working fluid.
# (He-Xe MW 83.8 was chosen to equal krypton's MW, so Kr and He-Xe share scalars.)
_MW_ARGON = 39.948


def gas_design_scalars(gas: str):
    """Return dict of compressor design scalars for the working fluid.
    We_des_lbps : design equivalent weight flow (= argon 0.581 * sqrt(MW/MW_argon))
    N_des_rpm   : design *physical* rotative speed at the 100% line
    U_tip_des_fps : design impeller tip speed (physical)
    """
    mw = {"argon": 39.948, "krypton": 83.80, "hexe_83p8": 83.8}[gas]
    f = (mw / _MW_ARGON) ** 0.5
    return dict(
        MW=mw,
        sqrtMW_ratio=f,
        We_des_lbps=COMP.We_des_lbps * f,           # argon 0.581 -> He-Xe/Kr ~0.841
        N_des_rpm=52200.0 / f,                       # argon 52200 -> He-Xe/Kr ~36000
        U_tip_des_fps=968.0 / f,                     # argon 968   -> He-Xe/Kr ~668
        PR_des=COMP.PR_des,                          # gas-independent
        eff_des_aero=COMP.eff_des,                   # rig (aerodynamic) design eff
    )


if __name__ == "__main__":
    # Self-check against Table II
    W, T1, P1 = 0.785, 540.0, 20.25
    print(f"theta={theta(T1):.4f}  delta={delta(P1):.4f}")
    print(f"W_corr = {w_corr(W, T1, P1):.4f} lb/s   (report 0.581)")
    print(f"N_corr = {n_corr(52200, T1):.0f} rpm     (report 51176)")
    print(f"U_corr = {n_corr(968, T1):.1f} ft/s      (report 949)")
