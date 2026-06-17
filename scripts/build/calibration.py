"""
Real-loop calibration transforms for the BRU compressor map.

Two effects the rig map does NOT contain, evidenced by the loop reports:

(1) HEAT-SOAK / "indicated" efficiency  (NASA TM X-67989, TM X-52826)
    The insulated, hot BRU soaks heat into the compressor, raising the MEASURED discharge
    temperature, so the temperature-rise-based ("indicated") efficiency reads low (loop
    69-75% vs rig 80-82%). Modeled as a fixed parasitic temperature rise dT_soak added to the
    aero temperature rise:
        dT_aero/T1 = (PR^((g-1)/g) - 1)/eta_aero            (actual aero rise, g=5/3 -> exp 0.4)
        eta_indicated = (PR^0.4 - 1) / ( dT_aero/T1 + dT_soak/T1 )
    dT_soak/T1 is a CALIBRATION constant; default 0.033 (~18 R at T1=540 R) puts peak indicated
    efficiency ~0.73, centered in the measured 69-75% band. Penalty is larger at low PR/low
    work (low-speed lines drop more), which is physically correct.

(2) HARDWARE / INSTRUMENTATION uncertainty BAND  (NASA TM X-52826)
    Evidence: He-Xe vs Kr indistinguishable, the three units indistinguishable; compressor PR
    ref/BRU 1.90 vs gas-loop 1.94. The apparent 74<->80% installation spread is OPERATING POINT
    (open-throttle vs design), already in the map. So the residual hardware+instrument band is
    SMALL: default +/- 1.5 efficiency points and +/- 2% PR. (nominal / optimistic / pessimistic)
"""
import numpy as np

GAMMA_EXP = 0.4            # (g-1)/g for g = 5/3
DT_SOAK_DEFAULT = 0.033   # dT_soak / T1  (calibrated to peak indicated eta ~0.73)
BAND_DEFF = 0.015         # +/- efficiency points (absolute) for hardware/instrument band
BAND_FPR = 0.02           # +/- fractional PR for the band


def indicated_efficiency(eta_aero, PR, dt_soak=DT_SOAK_DEFAULT):
    """Convert aerodynamic efficiency to heat-soak 'indicated' efficiency."""
    eta_aero = np.asarray(eta_aero, dtype=float)
    PR = np.asarray(PR, dtype=float)
    ideal = PR ** GAMMA_EXP - 1.0
    dt_aero = ideal / eta_aero
    return ideal / (dt_aero + dt_soak)


def apply_band(eff, PR, sign, d_eff=BAND_DEFF, f_PR=BAND_FPR):
    """sign=+1 optimistic, -1 pessimistic, 0 nominal. Returns (eff_b, PR_b)."""
    return np.asarray(eff) + sign * d_eff, np.asarray(PR) * (1.0 + sign * f_PR)


if __name__ == "__main__":
    # peak indicated efficiency at the design speed line
    for dt in (0.025, 0.033, 0.040):
        e = indicated_efficiency(0.819, 1.90, dt)   # peak aero eta at design PR-ish
        print(f"dT_soak/T1={dt:.3f} ({dt*540:.0f} R): peak indicated eta ~ {float(e):.3f}")
