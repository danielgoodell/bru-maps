"""
Compressor Reynolds-number correction for the BRU 4.25-in centrifugal compressor.

Source: NASA TN D-6640 (1972), "Reynolds Number Effect on Overall Performance of a
10.8-cm Sweptback-Bladed Centrifugal Compressor."

Definition (report): Re_U = rho1 * U_t3 * D_t3 / mu1
    rho1 = inlet total density, U_t3 = impeller exit tip speed, D_t3 = exit tip diameter,
    mu1  = inlet dynamic viscosity.   Design Re_U = 3.1e6.

Loss law (report eq.):   (1 - eta) / (1 - eta_ref) = [ (Re_U)_ref / Re_U ] ^ n
Reference (near design): Re_ref = 3.03e6, loss_ref = 1 - eta_max,ref = 0.187 (eta=0.813).
Piecewise exponents (report, three-part fit; n rises as Re falls):
    n = 0.06   for 1.21e6 < Re_U < 3.03e6
    n = 0.09   for 0.43e6 < Re_U < 1.21e6
    n = 0.20   for 0.34e6 < Re_U < 0.43e6
(single n = 0.1 approximates the whole tested range; above ~1e6 loss flattens, n -> 0.)

Self-check vs report: chaining the bands from 3.03e6 down to 0.34e6 gives loss 0.187 -> 0.227,
i.e. eta_max 0.813 -> 0.773 (a 4.0-point drop), matching the report's "1.5 + 2.5 points".

Secondary low-Re effects the report notes (applied as documented approximations):
  - peak PR falls 1.96 -> 1.90 (~3%) from Re 3.03e6 -> 0.34e6
  - max-flow and surge points shift to LOWER flow as Re drops
"""
import numpy as np
from corrected_params import P_STD_PSIA  # noqa: F401  (kept for unit context)

R_UNIV = 8.314462          # J/mol-K
D_T3_M = 0.108             # impeller exit tip diameter, m (4.25 in)

# loss-law anchors / bands  (upper Re of band, exponent within band)
RE_REF = 3.03e6
LOSS_REF = 0.187           # 1 - eta_max at Re_ref
_BANDS = [(3.03e6, 0.06), (1.21e6, 0.09), (0.43e6, 0.20), (0.34e6, 0.20)]

# Gas viscosity near 300 K for computing Re_U of a real loop condition.
# Sourced from the cycle project's NobleGasMixture backend (Tournier/El-Genk/Gallo virial EOS +
# transport correlations). CROSS-VALIDATED: argon mu=2.294e-5 reproduces the NASA TN D-6640 argon
# Re_U=3.07e6 at the 20-psia rig inlet (report design 3.1e6).
# mu is weakly T-dependent (~T^0.7); these are mu(300 K). For loop T in 289-322 K the spread is
# ~+/-4%; query ../cycle (GasCycle NobleGasMixture.viscosity) for exact values at other T.
GAS = {
    #            MW [kg/mol]   mu(300K) [Pa-s]   source
    "argon":     dict(MW=0.039948, mu=2.294e-5, flag="GasCycle-verified"),
    "krypton":   dict(MW=0.083800, mu=2.611e-5, flag="GasCycle"),
    "hexe_83p8": dict(MW=0.083800, mu=2.468e-5, flag="GasCycle"),
}


def density(P_pa, T_k, gas):
    return P_pa * GAS[gas]["MW"] / (R_UNIV * T_k)


def reynolds_U(P_pa, T_k, U_tip_mps, gas):
    """Compressor tip Reynolds number Re_U = rho1 U_t3 D_t3 / mu1."""
    rho = density(P_pa, T_k, gas)
    return rho * U_tip_mps * D_T3_M / GAS[gas]["mu"]


def loss_scale(Re_U):
    """S(Re) = (1-eta)/(1-eta_ref): multiplicative loss factor relative to Re_ref.
    Continuous, chained across the piecewise-n bands.  S(Re_ref)=1; S>1 below ref."""
    Re_U = np.asarray(Re_U, dtype=float)

    def one(re):
        if re >= RE_REF:                       # loss flattens at high Re
            return (RE_REF / re) ** 0.06
        S, top = 1.0, RE_REF
        for lo, n in [(1.21e6, 0.06), (0.43e6, 0.09), (0.34e6, 0.20), (0.0, 0.20)]:
            seg_lo = max(re, lo)
            if seg_lo < top:
                S *= (top / seg_lo) ** n
                top = seg_lo
            if re >= lo:
                break
        return S
    return np.vectorize(one)(Re_U)


def correct_efficiency(eta_map, Re_U):
    """Map efficiency (at Re_ref) -> efficiency at Re_U.  eta = 1 - (1-eta_map)*S(Re)."""
    return 1.0 - (1.0 - np.asarray(eta_map)) * loss_scale(Re_U)


def correct_PR(PR_map, Re_U):
    """Approximate PR knock-down at low Re, tied to the loss increase.
    Calibrated so peak PR 1.96 -> 1.90 (factor 0.969) at Re=0.34e6 (S~1.215)."""
    S = loss_scale(Re_U)
    # k chosen so (1 - k*(S-1)) = 0.969 when S = loss_scale(0.34e6)
    S_lo = float(loss_scale(0.34e6))
    k = (1 - 0.969) / (S_lo - 1)
    return np.asarray(PR_map) * (1.0 - k * (S - 1.0))


if __name__ == "__main__":
    # 1) argon reproduces the report Re at the base-map inlet (20.25 psia, 540 R, U=295 m/s)
    P = 20.25 * 6894.76
    Re = reynolds_U(P, 300.0, 295.0, "argon")
    print(f"argon base-map Re_U = {Re:.3e}  (report design ~3.1e6)")
    # 2) loss-law self-check: eta_max from 3.03e6 down to 0.34e6
    for re in [3.03e6, 1.21e6, 0.43e6, 0.34e6]:
        eta = 1 - LOSS_REF * float(loss_scale(re))
        print(f"  Re={re:.2e}  S={float(loss_scale(re)):.4f}  eta_max={eta:.4f}")
