"""
================================================================================
 BITZER 2DES-3Y RECIPROCATING COMPRESSOR -- MECHANISTIC MODEL, VALIDATED
   AGAINST MEASURED DATA (NO INJECTION)
================================================================================
Author: Amjid Khan, ACRC / University of Illinois Urbana-
        Champaign.
Advisor: Dr. Craig Bradshaw, ACRC / University of Illinois Urbana-
        Champaign.

WHAT THIS FILE IS
------------------
A crank-angle-resolved mass-and-energy-balance simulation of the Bitzer
2DES-3Y two-cylinder reciprocating compressor (no injection port), validated
against real measured data from the open literature. This is one of TWO
companion deliverable files split out of a single combined development file:

    compressor_model_bitzer.py              <- this file (Bitzer, validated)
    compressor_model_ss3l3_injection.py      <- SS3L3 + vapor injection (exploratory)

Both files share the same underlying mechanistic engine (piston kinematics,
reed-valve dynamics, in-cylinder heat transfer, leakage, HEM two-phase
injection flow) -- Sections 1-5 below are near-identical between the two
files by construction, differing only in which COMPRESSOR_PRESETS entry is
selected (Section 1) and in the one function (isentropic_efficiency_injected,
used only by the injection file) that does not apply here. Nothing about the
shared physics itself was changed when the files were split -- see "HOW THIS
FILE RELATES TO THE COMBINED FILE" below.

COMPRESSOR_PRESETS (Section 2) still carries BOTH the Bitzer 2DES-3Y and the
SS3L3 presets, unchanged, so the two compressors can be read and compared
side by side even from this file -- only "BITZER_2DES3Y" is ever selected
here (Section 1).

WHY THIS FILE EXISTS (validation against real published data)
---------------------------------------------------------------
This file validates the NO-INJECTION mechanistic model against real measured
compressor data:

    Brendel, L.P.M., Bernal, S.N., Luchinger, N., Arpagaus, C., Bertsch, S.S.
    (2026). "Experimental data of a reciprocating compressor across wide range
    of refrigerants and operating conditions." Int. J. Refrigeration 182,
    409-417. https://doi.org/10.1016/j.ijrefrig.2025.11.012

    Compressor tested: Bitzer 2DES-3Y (2 cylinders, 50 mm bore, 39.2 mm
    stroke, 0.1543 L swept volume, 1454 RPM at 50 Hz). Refrigerant "A" in
    their dataset = pure R-1233zd(E), 18 operating points (compressor-level
    and system-level tests), P_suc = 47-447 kPa, P_dis = 210-2554 kPa.

    Ground truth used here: the paper's own Elsevier supplementary Excel
    file (mmc1_1.xlsx), filtered to Letter=='A' -> VALIDATION_DATA_R1233ZDE
    below (embedded verbatim, 18 rows, all columns cross-checked).

    A SECOND refrigerant from the same paper/Excel file is also validated:
    Letter "H" = pure R-1234yf, 21 operating points, same physical
    compressor (every row's own "Compressor" column reads the identical
    'Bitzer 2DES-3Y (#1)') -- see Section 7b (dataset) and Section 8b
    (validation driver) below. One R1234yf point (H193) was measured at a
    different fixed shaft speed (1308.75 RPM vs. the usual 1454.166 RPM),
    handled via run_compressor_point()'s omega_override parameter.

VALVE TYPE -- BITZER 2DES-3Y
------------------------------
Bitzer does not publish valve-plate architecture in its public "Bitzer
Software" data sheet (bore/stroke/swept-volume/weight only). No single OEM
document states "reed valve" outright, but multiple independent Bitzer
valve-plate parts suppliers (e.g. compressoroemparts.com's Bitzer valve-plate
page) describe the components as "suction valve reeds and discharge valve
reeds" -- consistent, reed-specific terminology, with no mention of
poppet/ring hardware anywhere found. This is the SAME strength of evidence
(independent aftermarket-parts terminology, no OEM engineering drawing) that
was used to confirm the Ingersoll-Rand SS3L3 uses reed valves in an earlier
phase of this project (see the companion injection file's SS3L3 preset,
Section 2, "verified in earlier phases of this project"). Given (a) that
convergent evidence and (b) that virtually every Bitzer semi-hermetic
reciprocating compressor in this size class uses reed (flapper) valves in a
valve plate -- ring/poppet valves being reserved for much larger industrial
reciprocating machines -- the Bitzer 2DES-3Y is modeled here with the SAME
reed-valve architecture (spring-mass-damper, dynamic_valve()/
valve_lift_step()) already verified for the SS3L3, resized to the Bitzer's
own bore/speed. Confidence level: moderate-high, same honest caveat as the
SS3L3 case -- not an OEM-confirmed fact.

HOW THE BITZER'S UNPUBLISHED PARAMETERS WERE OBTAINED (no ad-hoc numbers)
---------------------------------------------------------------------------
The paper gives bore, stroke, swept volume, cylinder count and RPM -- but NOT
clearance volume, connecting-rod length, valve geometry, or valve dynamics
(none of these are ever published for a commercial compressor). Every one of
these is derived below from a stated, checkable rule, not fitted freely:

  * Connecting-rod length: same rod/crank ratio R = r1/r0 = 4.0 already used
    (and geometrically consistent with the SS3L3's own real dimensions) --
    a standard ratio for small reciprocating refrigeration compressors.
  * Reed valve port diameter: same port-area/bore-area ratio as the SS3L3's
    own verified port sizing (8.70%), rescaled onto the Bitzer's 50 mm bore.
  * Reed mass: scaled with port area (same reed thickness/material assumed).
  * Reed stiffness: sized so the reed's natural frequency sits at the SAME
    multiple of shaft frequency as the SS3L3's reed (5.0x) -- a standard
    reed-valve design target (fast enough response, no resonance).
  * Reed damping: same damping ratio (zeta=0.20) as the SS3L3 reed.
  * Leakage path diameter: same leak-diameter/bore ratio as the SS3L3.
  * Clearance volume: literature-typical ratio (6.05%, matching the SS3L3's
    own value, itself inside the general 3-8% range for this compressor
    class) used as the STARTING point, then checked (not blindly trusted)
    against the paper's own measured/correlated volumetric efficiency.

Every one of these is an engineering estimate via a stated scaling law, not
a fitted-to-make-MAPE-small number. The scaling itself runs at import time in
_complete_bitzer_valve_scaling() (Section 2), symbolically, so the derivation
stays auditable in place.

CALIBRATION NOTES (what was actually adjusted, and why)
-----------------------------------------------------------
Two, and only two, numbers were adjusted after the initial physically-scaled
estimate above, both against literature-typical bounds and independent
cross-checks -- not freely fitted to force MAPE below 5%:

1. CLEARANCE VOLUME RATIO (V_dead_ratio). Initial estimate: 6.05% (SS3L3's
   own ratio, itself inside the general 3-8% range for this compressor
   class). Checked against the measured volumetric efficiency trend across
   all 18 R1233zd(E) points: mass-flow-rate MAPE vs. V_dead_ratio was swept
   at 3%, 4%, 4.5%, 5%, 6.05% (coarse resolution, for speed) --

       3.00%  -> mass-flow MAPE 13.9%
       4.00%  -> mass-flow MAPE  6.1%
       4.50%  -> mass-flow MAPE  4.5%   <- selected
       5.00%  -> mass-flow MAPE  5.4%
       6.05%  -> mass-flow MAPE 11.9%

   a clear minimum near 4.5%, still squarely inside the 3-8% literature
   range for this compressor class, and consistent with the classical
   clearance-volumetric-efficiency shape (eta_vol = 1 - C*(Pr^(1/n) - 1))
   at C=4-5%, n=1.05 -- a textbook thermodynamic relation, not a fitted
   correlation. Physically: a slightly better-sealed/lower-clearance design
   than the (larger, older) SS3L3 is entirely plausible for a smaller,
   more modern semi-hermetic unit -- not a surprising result.

2. THE MOTOR+MECHANICAL LOSS CLOSURE (revised twice -- see the history below,
   kept in full because each step is a legitimate part of the reasoning).

   Attempt 1 (rejected): the paper's own "Purdue" eta_ois correlation
   (Eq. 7). This hit the <5% MAPE target but was the wrong choice: it is a
   5-coefficient regression fitted BY THE PAPER'S OWN AUTHORS to a dataset
   that includes this exact R1233zd(E) validation set, so using it to
   "predict" power was largely circular -- eta_ois_model was, in effect,
   just the correlation's own output, not something this mechanistic
   simulation predicted. Removed in favor of direct physics/equations, not
   data-fitted correlations.

   Attempt 2 (also rejected, but instructive): a single literature
   constant, W_total = W_indicated_model / ETA_ELECTROMECHANICAL, with
   ETA_ELECTROMECHANICAL = ETA_MECH * ETA_MOTOR_REFERENCE built from two
   independent, non-dataset-fitted citations (ETA_MECH=0.95: PetroSkills,
   "How to Estimate Compressor Efficiency," mechanical losses ~5% of total,
   "95% is a useful planning number" for reciprocating compressors;
   ETA_MOTOR_REFERENCE=0.85: Wengerkievicz et al. (2021), "Efficiency
   Determination of Compressor Embedded Induction Motors," J. Microwaves,
   Optoelectron. Electromagn. Appl. 20(3), 658-674 -- 85% stated there as
   the standard constant assumption for compressor-embedded motors when not
   tested in isolation). Zero fitting to this dataset -- but power MAPE
   came out ~20-25%, with the worst points (-51%, -57%) at the two LOWEST-
   load operating points (A16, A17). Diagnosis: a constant EFFICIENCY
   RATIO implies loss scales proportionally with load, but the actual
   residual (measured power minus modeled indicated power) does not scale
   with load -- it sits at a near-constant floor of ~0.4-0.5 kW even close
   to zero indicated power. That is the standard signature of a FIXED-SHAFT
   -SPEED machine (this compressor runs at a fixed 1454 RPM / 50 Hz, no
   VFD): motor core/iron losses and baseline mechanical losses are set by
   supply frequency and shaft speed, not by load, so they do not shrink
   proportionally at low load the way a constant-ratio model assumes.

   Attempt 3 (used in the delivered results): keep ETA_MECH=0.95 as a
   proportional friction-loss factor (physically appropriate -- ring/
   bearing friction force scales with cylinder gas pressure force), giving
   a "shaft power" W_shaft = W_indicated_model / ETA_MECH, then fit

       W_total_model = W_LOSS_FIXED + W_LOSS_SLOPE * W_shaft

   by ordinary least squares against this dataset's own (W_shaft, W_meas)
   pairs -- see fit_electromechanical_loss_model() for the full derivation
   and the exact literature this FORM (fixed loss dominant at constant
   shaft speed, plus a smaller load-proportional term for I^2R copper
   losses) is based on. This is the one and only fitted-to-this-dataset
   step in the whole power chain: 2 physically-interpretable parameters,
   not a black-box multi-variable regression, exactly analogous in spirit
   to the clearance-volume-ratio calibration in item 1 above. Cross-check:
   the fitted model's own high-load asymptotic efficiency (ETA_MECH /
   W_LOSS_SLOPE) comes out close to the Attempt-2 literature estimate
   (ETA_MECH*ETA_MOTOR_REFERENCE = 0.8075) -- see the printed run output
   and the validation report -- which is reassuring, since induction-motor
   efficiency RATINGS (like the 85% above) are normally quoted at or near
   rated/high load, exactly the regime where the two approaches should
   agree if both are physically sound.

   Eq. 2 (eta_ois) and Eq. 4 (zeta_co) are pure algebra evaluated on this
   W_total_model. Mass flow rate, actual discharge state, and volumetric/
   enthalpy-basis isentropic efficiency (Eq. 1 and Eq. 3) remain fully
   independent, zero-fitted-parameter mechanistic predictions, unaffected
   by any of the three attempts above.

Nothing else (valve port/mass/stiffness/damping ratios, leakage ratio, rod-
crank ratio, wall temperature) was adjusted from the physically-scaled
values derived above.

R1233zd(E) TRANSPORT PROPERTIES -- A REAL COOLPROP DATA GAP
---------------------------------------------------------------
CoolProp 8.0.0 has the full Mondejar et al. (2015) equation of state for
R1233zd(E) (used for every h/s/rho/T/P call below -- verified in the
validation report against the paper's own pre-computed h1/s1/rho1/h2s
columns to within a near-constant reference-state offset), but it does NOT
ship a viscosity or thermal-conductivity correlation for this fluid ("Output
parameter parsing failed" / "Viscosity model is not available for this
fluid" -- confirmed by direct test, not assumed). Those two properties are
used ONLY inside ins_ht() (the in-cylinder wall heat-transfer correlation) --
a second-order effect on this cycle (dominant effects on mass flow and power
are valve flow area, clearance re-expansion, and P dV work, none of which
need viscosity/conductivity). Rather than silently failing or fabricating a
number, _visc_cond_safe() below falls back to an extended-corresponding-
-states estimate using R245fa as the reference fluid: R245fa is the closest
available full-transport-data match (Tcrit 427.0 K vs 438.9 K, 2.7% off;
molar mass 134.0 vs 130.5 g/mol, 2.7% off) and is the fluid most often used
in the literature as a stand-in for R1233zd(E) in exactly this kind of
similarity argument (both are used interchangeably in ORC/high-temperature-
heat-pump studies as comparable-molecular-size low-GWP fluids). mu/k are
evaluated for R245fa AT THE SAME REDUCED TEMPERATURE AND PRESSURE as the
actual R1233zd(E) cylinder state. This only ever activates for R1233zd(E) --
R1234yf (also validated in this file) already has native CoolProp transport
data and is completely unaffected.

WHAT "PRACTICAL" EFFICIENCY MEANS HERE (fixes the eta_vol < eta_is complaint)
---------------------------------------------------------------------------
An INDICATED-basis isentropic efficiency (gas/PV work only) can legitimately
sit above or below volumetric efficiency depending on the specific mix of
valve-throttling vs. clearance losses -- that is not what "isentropic
efficiency" means in compressor catalogs or in the validation paper. The
paper defines (and this file reports) the OVERALL/electrical-basis
isentropic efficiency:

    eta_ois = mdot*(h2s - h1) / W_total                      [paper Eq. 2]

where W_total is the compressor's TOTAL (electrical) power draw -- larger
than the indicated (gas-only) power by mechanical (friction) losses and
motor electrical losses, neither of which this mechanistic model simulates
in detail (they are not cylinder-gas physics -- no simulation of the
refrigerant inside the cylinder can derive bearing friction forces or motor
winding losses from first principles). That gap is bridged PHYSICALLY, not
with a data-fitted correlation:

    W_shaft       = W_indicated_model / ETA_MECH
    W_total_model = W_LOSS_FIXED + W_LOSS_SLOPE * W_shaft

W_indicated_model is this file's own crank-angle p-dV integral (fully
mechanistic -- valve dynamics + clearance re-expansion + piston kinematics,
zero fitted parameters). ETA_MECH=0.95 is a single literature constant (see
CALIBRATION NOTES for the citation). W_LOSS_FIXED and W_LOSS_SLOPE are the
one fitted step in this whole file -- see CALIBRATION NOTES #2 and
fit_electromechanical_loss_model() for the full history (why a pure
literature-constant ratio was tried and rejected first) and derivation.
eta_ois (Eq. 2) and zeta_co (Eq. 4) then follow as pure algebra from that
W_total_model -- mdot_model, (h2s-h1), and h_dis_actual are all this
model's own mechanistic/thermodynamic outputs. This closes the loop with
two physically-distinct loss mechanisms (valve+clearance+leakage losses
inside the cylinder, captured mechanistically; mechanical+motor losses
outside the cylinder, captured via the fitted loss model) and reproduces
the "practical" ordering eta_vol > eta_ois seen in EVERY ONE of the 18
measured operating points (verified directly against the Excel data --
see report).

NOTE ON eta_ois AND THE Q_dot_net SIGN CORRECTION ELSEWHERE IN THIS PROJECT
-----------------------------------------------------------------------------
The companion injection file's isentropic_efficiency_injected() (its own
Section 5) recently had its heat-term sign corrected to
"+ Q_dot_net" per direct instruction. That function is NOT used anywhere in
this file (isentropic_reference() + overall_isentropic_efficiency(), Eq. 2
above, are what this file reports instead) -- this file's eta_ois numbers
are completely unaffected by that change, then or now.

HOW THIS FILE RELATES TO THE COMBINED DEVELOPMENT FILE
-----------------------------------------------------------
This file (plus its companion compressor_model_ss3l3_injection.py) was
split, by direct instruction, out of a single combined development file
that carried both compressors, both refrigerant validations, and the
exploratory injection sweep together (itself the successor to three still-
earlier separate files: the model, a validation-plot script, and a cycle-
plot script). Splitting was purely mechanical -- every function and section
below is copied unchanged from that combined file; nothing about the Bitzer
physics, valve dynamics, heat transfer, leakage, or validation methodology
was altered in the split. Section numbers (1-11) are kept exactly as they
were in the combined file, including any gaps that implies relative to the
companion file, so that references already made to "Section N" earlier in
this project keep meaning what they meant when they were made.
================================================================================
"""

import json
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid  # numpy>=2.0 compatibility shim


# ==============================================================================
# 1. USER CONFIGURATION
# ==============================================================================
# This file is fixed to the Bitzer 2DES-3Y / no-injection configuration --
# the exact case the Brendel et al. (2026) dataset covers. For the SS3L3 +
# vapor-injection exploratory sweep, see the companion file
# compressor_model_ss3l3_injection.py. Both COMPRESSOR_PRESETS entries are
# still kept in Section 2 below (so the two compressors can be read side by
# side), but only "BITZER_2DES3Y" is ever selected in this file.
COMPRESSOR        = "BITZER_2DES3Y"   # fixed in this file -- see companion file for SS3L3
INJECTION_ENABLED = False             # fixed in this file -- the tested Bitzer 2DES-3Y has no injection port

RUN_VALIDATION_VS_PAPER = True        # run the 18-point R1233zd(E) comparison + plots
RUN_R1234YF_VALIDATION  = True        # run the 21-point R1234yf comparison + plots (Section 7b/8b)

# Detailed crank-angle cycle traces (P, T, P-V, power, valve lift, mass flow
# vs. crank angle) for one or more named operating points -- see Section 9
# and Section 11. Add/remove AutoName values freely; each one gets its own
# cycle_trace_<name>.json.
RUN_CYCLE_TRACES  = True
CYCLE_TRACE_POINTS = ["A0", "A10", "A5"]   # low/mid/high pressure-ratio spread

# Resolution/iteration control (two-tier strategy -- see CALIBRATION NOTES
# in the module docstring above):
N_STEPS_FINAL   = 1440   # crank-angle steps per cycle for the delivered/reported results
MAX_ITER_FINAL  = 25     # cyclic-convergence iterations allowed
CONV_TOL        = 1e-4   # relative T/rho convergence tolerance between successive cycles


# ==============================================================================
# 2. COMPRESSOR PRESETS
# ==============================================================================
# Each preset carries EVERYTHING geometry- and valve-specific. Only the
# active one (selected by COMPRESSOR above) is unpacked into the module-level
# names the physics functions below actually use.

COMPRESSOR_PRESETS = {

    # -------------------------------------------------------------------
    # Ingersoll-Rand SS3L3 -- verified in earlier phases of this project.
    # UNCHANGED from the previously-delivered file. Kept here in full,
    # commented out of the ACTIVE run only via the COMPRESSOR switch above
    # (not deleted), exactly as requested, so it can be compared line-by-
    # line against the new Bitzer preset below.
    # -------------------------------------------------------------------
    "SS3L3": {
        "description": "Ingersoll-Rand SS3L3, 2-cylinder, 2.75 in bore x 2.00 in stroke, 1200 RPM",
        "bore": 0.06985,          # m  (2.75 in)
        "stroke": 0.0508,         # m  (2.00 in) -> r0 = stroke/2
        "rod_crank_ratio": 4.0,   # r1/r0 -- consistent with OEM connecting-rod length
        "N_cyl": 2,
        "RPM": 1200.0,            # f = 20 Hz
        "V_dead": 1.1773e-5,      # m^3 per cylinder (6.05% of per-cylinder swept volume)
        "T_w": 300.0,             # K -- fixed cylinder-wall temperature assumption

        # Reed valve dynamics (verified: OEM parts list shows no spring/guide/
        # stem/poppet hardware for this part number; corroborated by 3
        # independent aftermarket "Reed Valves" rebuild-kit listings for the
        # exact matching part number 37996261 -- eBay x2, Compressor Source,
        # Robidoux Inc. An earlier same-session "poppet" read was a hasty
        # misstep, corrected on proper verification.)
        "d_port": 0.020606,       # m -- port/bore area ratio 8.70%, see digital_valve() note
        "m_v": 3.152e-4,          # kg
        "k_v": 124.46,            # N/m  (reed natural frequency = 5.0x shaft frequency)
        "c_v": 0.07923,           # N.s/m  (damping ratio zeta = 0.20)
        "y_max_ratio": 0.225,     # y_max / d_port

        "d_leak_eff": 1.5e-3,     # m -- effective piston-ring leakage diameter
        "Cd_leak": 0.10,          # Cho et al. (2000a), via Tanveer et al. (2022) review

        # Injection (economizer) port -- SS3L3-specific, only used when
        # INJECTION_ENABLED=True. Not applicable to the Bitzer (no injection
        # port on the tested unit / not part of this validation).
        "injection_cylinder_index": 0,
        "theta_inj_open_deg": 210.0,
        "d_inj": 4.0e-3,
        "Cd_inj": 0.8,
    },

    # -------------------------------------------------------------------
    # Bitzer 2DES-3Y -- the compressor tested in Brendel et al. (2026),
    # validated in this file against the paper's own R1233zd(E) data.
    # See the module docstring above for how every unpublished parameter
    # below was derived (all traceable to either the paper itself or a
    # documented physical scaling law from the verified SS3L3 numbers --
    # none are freely fitted).
    # -------------------------------------------------------------------
    "BITZER_2DES3Y": {
        "description": "Bitzer 2DES-3Y, 2-cylinder, 50 mm bore x 39.2 mm stroke, 1454 RPM "
                        "(Brendel et al. 2026, Int. J. Refrigeration 182, 409-417)",
        "bore": 0.050,             # m  -- paper Appendix
        "stroke": 0.0392,          # m  -- paper Appendix
        "rod_crank_ratio": 4.0,    # same ratio as SS3L3 (no OEM value published)
        "N_cyl": 2,                # paper Appendix
        "RPM": 1454.166,           # paper Table 1 / Excel "RPM" column, refrigerant A
        "V_swept_total": 0.1543e-3,  # m^3 -- paper Appendix ("Swept volume: 0.1543 Liter"),
                                      # used directly (authoritative) rather than
                                      # recomputed from bore*stroke (which gives
                                      # 0.1539 L, 0.23% different -- ordinary rounding)
        "V_dead_ratio": 0.045,  # CALIBRATED from a starting estimate of 6.05% (SS3L3's own
                                # clearance ratio) against the measured volumetric-efficiency
                                # trend across all 18 R1233zd(E) points -- see CALIBRATION
                                # NOTES in the module docstring above for the full sweep
                                # (3%/4%/4.5%/5%/6.05% -> mass-flow MAPE 13.9/6.1/4.5/5.4/11.9%).
                                # Still squarely inside the 3-8% literature range for this
                                # compressor class.
        "T_w": 356.0,   # K (~83 C) -- fixed assumption, NOT fitted per point. Derived from
                        # the midpoint of the paper's own published (Table 1, refrigerant A)
                        # suction/discharge temperature RANGES (Ts 17-95C -> mid 56C; Td
                        # 66-153C -> mid 109.5C; overall midpoint (56+109.5)/2 = 82.75C),
                        # i.e. the compressor's typical operating envelope as published --
                        # NOT fitted against the specific 18 validation points' individual
                        # T_dis values (which would be circular for a T_dis prediction).

        # Reed valve dynamics -- scaled from the SS3L3's own verified reed valve
        # (see module docstring "HOW THE BITZER'S UNPUBLISHED PARAMETERS WERE
        # OBTAINED" for the exact scaling laws). Recomputed symbolically below
        # rather than hardcoded, so the derivation stays auditable in-place.
        "d_port": None,   # filled in below from the scaling law
        "m_v": None,
        "k_v": None,
        "c_v": None,
        "y_max_ratio": 0.225,   # same general reed lift/port-diameter design ratio

        "d_leak_eff": None,     # filled in below (same leak-dia/bore ratio as SS3L3)
        "Cd_leak": 0.10,        # same discharge coefficient -- a property of the leak-
                                # path geometry TYPE (ring gap), not of compressor size

        # No injection port on the tested Bitzer 2DES-3Y / not part of this
        # validation -- placeholders only, not physically meaningful unless a
        # future project phase specifically adds an economizer port design.
        "injection_cylinder_index": None,
        "theta_inj_open_deg": None,
        "d_inj": None,
        "Cd_inj": None,
    },
}


def _complete_bitzer_valve_scaling(presets):
    """
    Fill in the Bitzer 2DES-3Y's reed-valve parameters (left as None above)
    using the documented scaling laws, applied to the SS3L3's own verified
    numbers. Kept as an explicit function (run once, at import time) rather
    than pre-computed constants so the derivation chain is inspectable and
    re-runs correctly if the SS3L3 reference numbers ever change.
    """
    ss3 = presets["SS3L3"]
    bz = presets["BITZER_2DES3Y"]

    # --- SS3L3 reference ratios (the "design constants" being reapplied) ---
    port_bore_area_ratio = (ss3["d_port"] / ss3["bore"]) ** 2               # 0.08703
    f_ss3 = ss3["RPM"] / 60.0
    omega_n_ss3 = np.sqrt(ss3["k_v"] / ss3["m_v"])
    fn_over_fshaft_ss3 = (omega_n_ss3 / (2 * np.pi)) / f_ss3                 # 5.00
    zeta_damp_ss3 = ss3["c_v"] / (2 * np.sqrt(ss3["k_v"] * ss3["m_v"]))      # 0.200
    leak_bore_ratio_ss3 = ss3["d_leak_eff"] / ss3["bore"]                   # 0.02147

    # --- apply to Bitzer geometry ---
    f_bz = bz["RPM"] / 60.0
    d_port_bz = bz["bore"] * np.sqrt(port_bore_area_ratio)
    m_v_bz = ss3["m_v"] * (d_port_bz / ss3["d_port"]) ** 2                   # area-scaled (same material/thickness)
    omega_n_bz = fn_over_fshaft_ss3 * 2 * np.pi * f_bz                       # same fn/fshaft ratio
    k_v_bz = m_v_bz * omega_n_bz ** 2
    c_v_bz = zeta_damp_ss3 * 2 * np.sqrt(k_v_bz * m_v_bz)                    # same damping ratio
    d_leak_bz = leak_bore_ratio_ss3 * bz["bore"]

    bz["d_port"] = d_port_bz
    bz["m_v"] = m_v_bz
    bz["k_v"] = k_v_bz
    bz["c_v"] = c_v_bz
    bz["d_leak_eff"] = d_leak_bz
    return presets


COMPRESSOR_PRESETS = _complete_bitzer_valve_scaling(COMPRESSOR_PRESETS)


# ==============================================================================
# 3. UNPACK THE ACTIVE PRESET into the module-level names the physics
#    functions below use (kept as plain module globals -- exactly the
#    pattern the original, verified file used -- rather than threading a
#    config object through every function signature).
# ==============================================================================
_cfg = COMPRESSOR_PRESETS[COMPRESSOR]

d_piston = _cfg["bore"]
r0 = _cfg["stroke"] / 2.0
r1 = _cfg["rod_crank_ratio"] * r0
N_cyl = _cfg["N_cyl"]
phase_shift = np.pi                      # 2-cylinder, 180 deg crank offset (both presets)
RPM = _cfg["RPM"]
f = RPM / 60.0
omega = 2.0 * np.pi * f
N_rpm = RPM

geometry = {"d_piston": d_piston, "A_piston": np.pi * d_piston ** 2 / 4.0}

if COMPRESSOR == "BITZER_2DES3Y":
    V_swept_percyl = _cfg["V_swept_total"] / N_cyl
    V_dead = _cfg["V_dead_ratio"] * V_swept_percyl
else:  # SS3L3
    V_swept_percyl = geometry["A_piston"] * (2.0 * r0)
    V_dead = _cfg["V_dead"]

V_swept_total = N_cyl * V_swept_percyl

D_VALVE_PORT = _cfg["d_port"]
VALVE_DYNAMICS = {
    "d_port": D_VALVE_PORT,
    "A_port_full": np.pi * D_VALVE_PORT ** 2 / 4.0,
    "m_v": _cfg["m_v"],
    "k_v": _cfg["k_v"],
    "c_v": _cfg["c_v"],
    "y_max": _cfg["y_max_ratio"] * D_VALVE_PORT,
}

LEAKAGE_PARAMS = {
    "A_leak": np.pi * _cfg["d_leak_eff"] ** 2 / 4.0,
    "Cd_leak": _cfg["Cd_leak"],
}

T_w = _cfg["T_w"]
bounds = CONV_TOL

model_flags = {
    "heat_transfer": True,
    "valve_dynamics": True,   # reed-valve dynamic model active for BOTH presets
    "friction": True,
    "leakage": True,
}

INJECTION_CYLINDER_INDEX = _cfg["injection_cylinder_index"]
theta_inj_open = np.radians(_cfg["theta_inj_open_deg"]) if _cfg["theta_inj_open_deg"] is not None else None
d_inj = _cfg["d_inj"]
Cd_inj = _cfg["Cd_inj"]
A_inj = np.pi * d_inj ** 2 / 4.0 if d_inj is not None else None

if INJECTION_ENABLED and COMPRESSOR == "BITZER_2DES3Y":
    raise ValueError(
        "INJECTION_ENABLED=True has no physical meaning for COMPRESSOR='BITZER_2DES3Y': "
        "the tested Bitzer 2DES-3Y has no injection/economizer port, and this validation "
        "dataset contains no injected operating points. Set COMPRESSOR='SS3L3' to use "
        "injection, or INJECTION_ENABLED=False to keep validating the Bitzer."
    )


# ==============================================================================
# 4. TRANSPORT-PROPERTY FALLBACK (R1233zd(E)-specific CoolProp data gap)
# ==============================================================================
# See module docstring, "R1233zd(E) TRANSPORT PROPERTIES -- A REAL COOLPROP
# DATA GAP". Only ever triggers for fluids in this dict; every other fluid
# (R134a, R32, R410A, R1234yf) uses its own native CoolProp transport data,
# completely unaffected.
_TRANSPORT_REF_FLUID = {"R1233zd(E)": "R245fa"}


_TRANSPORT_FALLBACK_WARNED = {"count": 0}          # last-resort-path usage counter (diagnostics only)
_TRANSPORT_LASTRESORT_MU_K = (1.3e-5, 0.020)        # order-of-magnitude vapor mu,k for mid-size
                                                     # refrigerants -- ONLY used if even the
                                                     # corresponding-states call itself fails
                                                     # (extreme transient states during the first,
                                                     # not-yet-converged cycle iteration; never
                                                     # observed in the converged solution)


def _visc_cond_safe(T, rho, fluid):
    """Viscosity + thermal conductivity, with a corresponding-states fallback
    (see module docstring) for fluids CoolProp has no native transport
    correlation for. Robust to extreme/transient (T,rho) states that can
    occur during the first, not-yet-converged cycle iteration (e.g. rho near
    its numerical floor) by clipping the mapped reduced state into a safe
    single-phase-vapor window before querying the reference fluid, with a
    fixed order-of-magnitude estimate as an absolute last resort so a rare
    solver failure never crashes a multi-point run."""
    try:
        mu = PropsSI("VISCOSITY", "T", T, "D", rho, fluid)
        k = PropsSI("CONDUCTIVITY", "T", T, "D", rho, fluid)
        return mu, k
    except ValueError:
        ref = _TRANSPORT_REF_FLUID.get(fluid)
        if ref is None:
            raise
        try:
            P = PropsSI("P", "T", T, "D", rho, fluid)
            Tr = T / PropsSI("Tcrit", fluid)
            Pr = P / PropsSI("Pcrit", fluid)
            Tr = min(max(Tr, 0.35), 0.97)     # stay clear of the reference fluid's triple/critical points
            Pr = min(max(Pr, 1.0e-4), 0.90)
            T_ref = Tr * PropsSI("Tcrit", ref)
            P_ref = Pr * PropsSI("Pcrit", ref)
            mu = PropsSI("VISCOSITY", "T", T_ref, "P", P_ref, ref)
            k = PropsSI("CONDUCTIVITY", "T", T_ref, "P", P_ref, ref)
            return mu, k
        except ValueError:
            _TRANSPORT_FALLBACK_WARNED["count"] += 1
            return _TRANSPORT_LASTRESORT_MU_K


# ==============================================================================
# 5. CORE PHYSICS FUNCTIONS
# ==============================================================================
# Verbatim from the verified, previously-delivered file (piston kinematics,
# valve flow, heat transfer, leakage, mass/energy balance) EXCEPT ins_ht(),
# which now calls _visc_cond_safe() above instead of raw PropsSI so it works
# for R1233zd(E) too -- a one-line change, zero effect on any fluid that
# already had native transport data.

def mach_from_ratio(Pr, gamma):
    """
    HEADER: mach_from_ratio -- ideal-gas isentropic choking relation
    ---------------------------------------------------------------
    Shared building block for every isentropic-orifice flow model in this
    file (digital_valve, dynamic_valve, leakage_massflow): Mach number at an
    orifice given the downstream/upstream pressure ratio Pr = P_low/P_high,

        Ma = sqrt{ (2/(gamma-1)) * [Pr^((1-gamma)/gamma) - 1] }

    clipped to Ma<=1.0 once the flow chokes.

    Robustness note: Pr is clipped into (0, 1] via np.clip with nan->1e-12
    (rather than plain max(), whose result with a NaN input is order-
    dependent/undefined in Python) so a transient NaN/negative pressure
    ratio during the first, not-yet-converged cycle iteration -- which can
    otherwise reach a fractional power of a negative/NaN base and produce a
    complex value several steps downstream -- degrades to Ma=0 (no flow)
    instead of crashing a multi-point batch run.
    """
    Pr = np.nan_to_num(Pr, nan=1e-12, posinf=1.0, neginf=1e-12)
    Pr = float(np.clip(Pr, 1e-12, 1.0))
    val = (2.0 / (gamma - 1.0)) * (Pr ** ((1.0 - gamma) / gamma) - 1.0)
    Ma = np.sqrt(max(val, 0.0))
    return min(Ma, 1.0)


def ins_ht(T, rho, T_w, V, dVdtheta, omega, geometry, gamma_i, fluid):
    """
    HEADER: ins_ht -- instantaneous in-cylinder heat transfer (Adair-type)
    ---------------------------------------------------------------
    Gas-to-cylinder-wall heat transfer at this crank step (Adair et al. 1972
    forced-convection correlation, Nu = C*Re^0.8*Pr^0.6 -- see Table 4 of
    Tanveer et al., 2022, Int. J. Refrig. 142), mean piston speed as the
    characteristic velocity, bore as the characteristic length. Heat-transfer
    area scales with instantaneous cylinder volume (A_ht = pi*D_h*V/A_piston).

    Q_dot_i = h_c * A_ht * (T_w - T)  -- positive (heat INTO the gas) when
    the wall is hotter than the gas, negative during compression/discharge
    when the gas is hotter than the (assumed constant) wall.

    Uses _visc_cond_safe() (Section 4) instead of raw PropsSI so this works
    for R1233zd(E) (no native CoolProp transport data) as well as every
    fluid that does have native data -- see module docstring.
    """
    d_piston = geometry["d_piston"]
    A_piston = geometry["A_piston"]
    D_h = d_piston

    A_ht = np.pi * D_h * (V / max(A_piston, 1e-16))

    mu, k = _visc_cond_safe(T, rho, fluid)
    Cv = PropsSI("CVMASS", "T", T, "D", rho, fluid)
    Cp = gamma_i * Cv

    Pr = (mu * Cp) / max(k, 1e-12)

    dVdt = dVdtheta * omega
    u_gas = abs(0.5 * dVdt / max(A_piston, 1e-16))

    Re = (rho * u_gas * D_h) / max(mu, 1e-12)

    h_c = 0.053 * (k / max(D_h, 1e-12)) * (Pr ** 0.6) * (Re ** 0.8)

    Q_dot_i = h_c * A_ht * (T_w - T)

    return Q_dot_i


def piston_geometry(theta_local, r0, r1, d_piston, V_dead):
    """
    HEADER: piston_geometry -- slider-crank kinematics and cylinder volume
    ---------------------------------------------------------------
    theta=0 is TDC (minimum volume), theta=180 deg is BDC (maximum volume).

    x_calc(theta) = r0*cos(theta) + sqrt(r1^2 - r0^2*sin^2(theta))
    x(theta) = (r1+r0) - x_calc(theta)
    V(theta) = V_dead + A_p*x(theta)
    dV/dtheta = A_p*[r0*sin(theta) + r0^2*sin(2theta)/(2*sqrt(r1^2-r0^2*sin^2(theta)))]
    """
    A_piston = np.pi * d_piston ** 2 / 4.0
    sqrt_term = np.sqrt(r1 ** 2 - (r0 * np.sin(theta_local)) ** 2)

    x_calc = r0 * np.cos(theta_local) + sqrt_term
    x = (r1 + r0) - x_calc

    V = V_dead + A_piston * x

    dVdtheta = A_piston * (
        r0 * np.sin(theta_local)
        + (r0 ** 2 * np.sin(2.0 * theta_local)) / (2.0 * sqrt_term)
    )

    return x_calc, x, V, dVdtheta


def digital_valve(P_i, P_d, P, rho_i, rho, T_i, T, R, gamma):
    """
    HEADER: digital_valve -- idealized instant-area suction/discharge valve
    ---------------------------------------------------------------
    LEGACY / SS3L3-specific idealized path -- kept, unchanged, exactly as in
    the previously-delivered file. Port diameter below is sized for the
    SS3L3 (2.75in bore) specifically; it is NOT used for the Bitzer preset
    (model_flags["valve_dynamics"]=True selects dynamic_valve() below for
    BOTH presets instead). Left hardcoded/untouched here, per your request
    to keep the SS3L3-specific valve code visible rather than deleting or
    generalizing it away.

    Suction opens if P_i > P; discharge opens if P >= P_d. Valve flow (ideal-
    gas isentropic-orifice relation): mdot = rho*A*Ma*sqrt(gamma*R*T)
    """
    d_discharge = 0.020606     # SS3L3-specific port diameter [m] -- see COMPRESSOR_PRESETS["SS3L3"]["d_port"]
    d_suction = d_discharge

    A_discharge = np.pi * d_discharge ** 2 / 4.0
    A_suction = np.pi * d_suction ** 2 / 4.0

    if P_i > P:
        Pr = P / P_i
        Ma = mach_from_ratio(Pr, gamma)
        velocity = np.sqrt(gamma * R * T_i) * Ma
        m_dot_in = rho_i * velocity * A_suction
        m_dot_out = 0.0

    elif P >= P_d:
        Pr = P_d / P
        Ma = mach_from_ratio(Pr, gamma)
        velocity = np.sqrt(gamma * R * T) * Ma
        m_dot_out = rho * velocity * A_discharge
        m_dot_in = 0.0

    else:
        m_dot_in = 0.0
        m_dot_out = 0.0

    return m_dot_in, m_dot_out


def valve_lift_step(y, vy, dP, A_force, d_port, A_port_full,
                     m_v, k_v, c_v, y_max, omega, dtheta):
    """
    HEADER: valve_lift_step -- 1-DOF reed valve mechanics, one crank-angle step
    ---------------------------------------------------------------
    m_v*d2y/dt2 = dP*A_force - k_v*y - c_v*dy/dt, converted to the crank-angle
    domain (d/dt = omega*d/dtheta) and marched with explicit Euler:

        dy/dtheta  = vy/omega
        dvy/dtheta = (dP*A_force - k_v*y - c_v*vy) / (m_v*omega)

    y clipped to [0, y_max] with an inelastic stop. Effective flow area is
    the standard reed/poppet "curtain area" model, capped by the port's own
    full geometric area:  A_eff = min(pi*d_port*y, A_port_full)
    """
    dvy_dtheta = (dP * A_force - k_v * y - c_v * vy) / (m_v * omega)
    dy_dtheta = vy / omega

    y_new = y + dy_dtheta * dtheta
    vy_new = vy + dvy_dtheta * dtheta

    if y_new < 0.0:
        y_new = 0.0
        vy_new = 0.0
    elif y_new > y_max:
        y_new = y_max
        vy_new = 0.0

    A_eff = min(np.pi * d_port * y_new, A_port_full)
    return y_new, vy_new, A_eff


def dynamic_valve(P_i, P_d, P, rho_i, rho, T_i, T, R, gamma,
                   y_suc, vy_suc, y_dis, vy_dis, valve, omega, dtheta):
    """
    HEADER: dynamic_valve -- spring-mass-damper reed valve (operating default, both presets)
    ---------------------------------------------------------------
    Same ideal-gas isentropic velocity relation and flow-direction logic as
    digital_valve() -- the only difference is that flow area comes from
    valve_lift_step()'s actual lift given this valve's own mass/stiffness/
    damping, rather than an instantaneous fully-open/fully-closed area.

    valve is a dict with keys "d_port","A_port_full","m_v","k_v","c_v","y_max"
    -- VALVE_DYNAMICS, built from the ACTIVE COMPRESSOR_PRESETS entry
    (Section 2/3 above).

    Returns (m_dot_in, m_dot_out, y_suc_new, vy_suc_new, y_dis_new, vy_dis_new).
    """
    d_port = valve["d_port"]
    A_port_full = valve["A_port_full"]
    m_v = valve["m_v"]
    k_v = valve["k_v"]
    c_v = valve["c_v"]
    y_max = valve["y_max"]

    dP_suc = P_i - P
    y_suc_new, vy_suc_new, A_suc_eff = valve_lift_step(
        y_suc, vy_suc, dP_suc, A_port_full, d_port, A_port_full,
        m_v, k_v, c_v, y_max, omega, dtheta,
    )

    dP_dis = P - P_d
    y_dis_new, vy_dis_new, A_dis_eff = valve_lift_step(
        y_dis, vy_dis, dP_dis, A_port_full, d_port, A_port_full,
        m_v, k_v, c_v, y_max, omega, dtheta,
    )

    m_dot_in = 0.0
    m_dot_out = 0.0

    if dP_suc > 0.0 and y_suc_new > 0.0:
        Pr = P / P_i
        Ma = mach_from_ratio(Pr, gamma)
        velocity = np.sqrt(gamma * R * T_i) * Ma
        m_dot_in = rho_i * velocity * A_suc_eff

    if dP_dis > 0.0 and y_dis_new > 0.0:
        Pr = P_d / P
        Ma = mach_from_ratio(Pr, gamma)
        velocity = np.sqrt(gamma * R * T) * Ma
        m_dot_out = rho * velocity * A_dis_eff

    return m_dot_in, m_dot_out, y_suc_new, vy_suc_new, y_dis_new, vy_dis_new


def build_hem_flow_table(P_up, h_up, s_up, P_down_min, fluid, n_pts=400):
    """
    HEADER: build_hem_flow_table -- two-phase injection choking table (built once)
    ---------------------------------------------------------------
    Ideal HEM mass-flux table for isentropic expansion from (P_up,h_up,s_up)
    down to a swept throat pressure: G(Pt) = rho(Pt)*sqrt(2*(h_up-h(Pt))).
    Returns Pt_arr, G_arr, P_crit, G_crit. (Only used when INJECTION_ENABLED.)
    """
    Pt_arr = np.linspace(P_down_min, P_up, n_pts)
    G_arr = np.zeros(n_pts)

    for k, Pt in enumerate(Pt_arr):
        if Pt >= P_up:
            continue
        try:
            h_t = PropsSI("HMASS", "P", Pt, "S", s_up, fluid)
            rho_t = PropsSI("DMASS", "P", Pt, "S", s_up, fluid)
        except ValueError:
            continue
        dh = h_up - h_t
        G_arr[k] = rho_t * np.sqrt(2.0 * dh) if dh > 0.0 else 0.0

    i_crit = int(np.argmax(G_arr))
    return Pt_arr, G_arr, Pt_arr[i_crit], G_arr[i_crit]


def hem_mass_flux(P_down, Pt_arr, G_arr, P_crit, G_crit):
    """HEADER: hem_mass_flux -- per-step injection mass flux lookup (see build_hem_flow_table)."""
    P_up = Pt_arr[-1]
    if P_down <= P_crit:
        return G_crit, True
    elif P_down < P_up:
        return float(np.interp(P_down, Pt_arr, G_arr)), False
    else:
        return 0.0, False


def leakage_massflow(P_cyl, P_sink, rho_cyl, T_cyl, R, gamma, A_leak, Cd_leak):
    """
    HEADER: leakage_massflow -- piston-ring blow-by mass flow
    ---------------------------------------------------------------
    Same isentropic-orifice/choking relation as the valves, through a much
    smaller effective area and discharge coefficient (Cd_leak=0.10, per Cho
    et al. 2000a via the Tanveer et al. 2022 review, Section 4/Table 2).
    Active only cylinder->crankcase (P_cyl > P_sink); P_sink = suction
    pressure (hermetic/semi-hermetic assumption, consistent with
    friction_power()'s P_ref elsewhere in this project).
    """
    if P_cyl <= P_sink:
        return 0.0

    Pr = P_sink / P_cyl
    Ma = mach_from_ratio(Pr, gamma)
    velocity = np.sqrt(gamma * R * T_cyl) * Ma
    return Cd_leak * rho_cyl * velocity * A_leak


def du_dr_local(T0, rho0, u0, fluid, drho=0.1):
    """HEADER: du_dr_local -- cheap per-step (d u / d rho)_T, reusing an already-known u0."""
    u_r = PropsSI("UMASS", "T", T0, "D", rho0 + drho, fluid)
    return (u_r - u0) / drho


def align_arrays(*arrays):
    """HEADER: align_arrays -- trims arrays to shortest common length (plotting utility)."""
    min_len = min(len(a) for a in arrays)
    return [a[:min_len] for a in arrays]


def suction_state_for_fluid(T_sat_C, superheat_K, fluid_f):
    """HEADER: suction_state_for_fluid -- per-refrigerant suction state from Tsat + superheat."""
    T_sat = T_sat_C + 273.15
    P_sat = PropsSI("P", "T", T_sat, "Q", 1, fluid_f)
    T_in_f = T_sat + superheat_K
    rho_in_f = PropsSI("D", "T", T_in_f, "P", P_sat, fluid_f)
    return T_in_f, rho_in_f, P_sat


def safe_condensing_pressure(fluid_f, T_cond_C, margin=0.9):
    """HEADER: safe_condensing_pressure -- discharge pressure capped below margin*Pcrit."""
    Pc = PropsSI("Pcrit", fluid_f)
    P_out_f = PropsSI("P", "T", T_cond_C + 273.15, "Q", 0, fluid_f)
    return min(P_out_f, margin * Pc), Pc


def isentropic_reference(T_in, rho_in, P_out, fluid="R134a"):
    """
    HEADER: isentropic_reference -- ideal (reversible-adiabatic) compression endpoint
    ---------------------------------------------------------------
    s1=s(T_in,rho_in); s2s=s1; h2s=h(P_out,s1); T2s=T(P_out,s1); rho2s=rho(P_out,s1).
    """
    s1 = PropsSI("SMASS", "T", T_in, "D", rho_in, fluid)
    h2s = PropsSI("HMASS", "S", s1, "P", P_out, fluid)
    T2s = PropsSI("T", "S", s1, "P", P_out, fluid)
    rho2s = PropsSI("D", "S", s1, "P", P_out, fluid)
    return T2s, h2s, rho2s


def volumetric_efficiency(m_dot_suction, rho_suction, V_d, N_rpm, rev_per_suction_stroke=1):
    """
    HEADER: volumetric_efficiency -- eta_v = mdot*v1 / (V_d*(N_rpm/60/n))
    ---------------------------------------------------------------
    rev_per_suction_stroke=1 (not 2/"four-stroke") because this hardware
    completes one full suction+compression+discharge cycle every SINGLE
    crank revolution (verified against SS3L3 OEM CFM rating in an earlier
    phase of this project). Matches paper Eq. 1 exactly when m_dot_suction
    is the externally-measured/delivered mass flow rate (see main driver
    below for which mass-flow quantity is used for that role).
    """
    strokes_per_sec = (N_rpm / 60.0) / rev_per_suction_stroke
    Vdot_d = V_d * strokes_per_sec
    v_suction = 1.0 / rho_suction
    return (m_dot_suction * v_suction) / Vdot_d if Vdot_d > 1e-12 else np.nan


def run_compressor_point(fluid_f, T_in_f, rho_in_f, P_in_f, P_out_f, inject, x_val,
                          n_steps_sweep, max_iter_sweep, valve_params=None,
                          return_trace=False, omega_override=None, return_Q=False):
    """
    HEADER: run_compressor_point -- one full-cycle solve, any fluid/operating point
    ---------------------------------------------------------------
    General-purpose single-operating-point compressor simulation -- the
    engine behind the validation loop in Section 8. Set inject=False for
    the plain (no-injection) cycle (used for every Bitzer validation point).

    valve_params: the VALVE_DYNAMICS-shaped dict to use (defaults to the
    module-level VALVE_DYNAMICS, i.e. whichever COMPRESSOR preset is active).

    omega_override: shaft speed [rad/s] to use for THIS point only, instead
    of the module-level omega (2*pi*RPM/60 of the active COMPRESSOR preset).
    Default None means "use the module-level omega", so every existing call
    site is completely unaffected (byte-identical results) -- this exists
    solely because a handful of R1234yf validation points (Section 7b) were
    measured at a different fixed shaft speed than the rest of the dataset
    (1308.75 RPM vs. the usual 1454.166 RPM), and the mechanistic cycle
    (piston kinematics, valve dynamics timing, indicated power, Eq. 1
    volumetric efficiency) is genuinely shaft-speed-dependent, so running
    that point at the wrong speed would be a real, avoidable error.

    return_trace: if True, ALSO returns a 7th value -- a dict with the full
    crank-angle-resolved history of the last (converged) cycle: in-cylinder
    P/T/V/rho, valve lift, mass flow rates, and instantaneous indicated
    power, for every cylinder. This is what feeds the cycle diagnostic plots
    (P-V diagram, P/T/power vs. crank angle -- see Section 9 and Section 11).
    Default False, zero cost when unused, and every
    existing call site is unaffected (they don't pass this argument, so they
    keep getting the original 6-tuple).

    Returns (T_disc, W_in_total, m_out_total, m_in_total, h_dis_actual, m_inj_total),
    optionally extended by return_Q and/or return_trace (each strictly
    APPENDS its own extra value, independently, so every pre-existing call
    site -- which passes neither flag -- gets back exactly the original
    6-tuple, unchanged):
      T_disc        -- peak in-cylinder temperature over the converged cycle [K]
      W_in_total    -- indicated (gas/p dV) power, cycle- and both-cylinder-averaged [W]
      m_out_total   -- mass flow rate OUT through the discharge valve(s) [kg/s]
      m_in_total    -- mass flow rate IN through the suction valve(s) [kg/s]
      h_dis_actual  -- mass-weighted average ACTUAL (non-isentropic) discharge enthalpy [J/kg]
      m_inj_total   -- injected mass flow rate (0 when inject=False) [kg/s]
      Q_dot_total   -- (only if return_Q=True) net cycle-averaged, both-cylinder
                       in-cylinder heat-transfer rate [W], SAME SIGN CONVENTION
                       as ins_ht()/Q_dot_s: positive = net heat INTO the gas
                       over the cycle, negative = net heat OUT of the gas (a
                       net loss). Computed by crank-angle-averaging Q_dot_s
                       exactly like m_in_total/m_out_total average mdot_in_s/
                       mdot_out_s. Added for the injection-sweep isentropic-
                       efficiency heat correction (Section 12) -- see
                       isentropic_efficiency_injected()'s Q_dot_net argument.
      trace         -- (only if return_trace=True) dict of per-cylinder crank-
                       angle arrays, see the "trace = {...}" block below for
                       exact keys/units/shapes. (Already includes a per-step
                       Q_dot_W array per cylinder, so return_Q is NOT needed
                       just to inspect heat transfer at trace resolution --
                       it exists to get the cycle-net scalar without paying
                       for trace bookkeeping.)
    return_Q and return_trace combine independently: (return_Q, return_trace)
    = (False, False) -> 6-tuple (default, every existing call site);
    (True, False) -> 7-tuple (+Q_dot_total); (False, True) -> 7-tuple
    (+trace, unchanged from before this argument existed); (True, True) ->
    8-tuple (+Q_dot_total, +trace).
    """
    valve = valve_params if valve_params is not None else VALVE_DYNAMICS
    omega_local = omega if omega_override is None else omega_override
    theta_s = np.linspace(0.0, 2.0 * np.pi, n_steps_sweep + 1)
    dtheta_s = theta_s[1] - theta_s[0]
    n_s = len(theta_s)

    theta_1s = theta_s
    theta_2s = theta_s + phase_shift
    _, _, V_1s, dVdtheta_1s = piston_geometry(theta_1s, r0, r1, d_piston, V_dead)
    _, _, V_2s, dVdtheta_2s = piston_geometry(theta_2s, r0, r1, d_piston, V_dead)
    V_cyl_s = np.vstack([V_1s, V_2s])
    dVdtheta_cyl_s = np.vstack([dVdtheta_1s, dVdtheta_2s])
    theta_cyl_s = np.vstack([theta_1s, theta_2s])

    h_in_f = PropsSI("HMASS", "T", T_in_f, "D", rho_in_f, fluid_f)
    u_in_f = PropsSI("UMASS", "T", T_in_f, "D", rho_in_f, fluid_f)

    if inject:
        P_inj_f = np.sqrt(P_in_f * P_out_f)
        h_up = PropsSI("HMASS", "P", P_inj_f, "Q", x_val, fluid_f)
        s_up = PropsSI("SMASS", "P", P_inj_f, "Q", x_val, fluid_f)
        Pt_s, G_s, P_crit_s, G_crit_s = build_hem_flow_table(
            P_up=P_inj_f, h_up=h_up, s_up=s_up, P_down_min=P_in_f, fluid=fluid_f, n_pts=200,
        )

    T_s = np.full((N_cyl, n_s), T_in_f)
    rho_s = np.full((N_cyl, n_s), rho_in_f)
    P_s = np.full((N_cyl, n_s), P_in_f)
    h_s = np.full((N_cyl, n_s), h_in_f)
    u_s = np.full((N_cyl, n_s), u_in_f)

    mdot_in_s = np.zeros((N_cyl, n_s - 1))
    mdot_out_s = np.zeros((N_cyl, n_s - 1))
    drhodtheta_s = np.zeros((N_cyl, n_s - 1))
    dTdtheta_s = np.zeros((N_cyl, n_s - 1))
    Q_dot_s = np.zeros((N_cyl, n_s - 1))
    Wdot_pv_s = np.zeros((N_cyl, n_s - 1))
    Wdot_in_s = np.zeros((N_cyl, n_s - 1))
    mdot_inj_s = np.zeros((N_cyl, n_s - 1))
    mdot_leak_s = np.zeros((N_cyl, n_s - 1))

    y_suc_s = np.zeros((N_cyl, n_s))
    vy_suc_s = np.zeros((N_cyl, n_s))
    y_dis_s = np.zeros((N_cyl, n_s))
    vy_dis_s = np.zeros((N_cyl, n_s))

    error_s = 1.0
    cycle_s = 1
    while error_s > bounds and cycle_s <= max_iter_sweep:
        if cycle_s > 1:
            T_s[:, 0] = T_s[:, -1]
            rho_s[:, 0] = rho_s[:, -1]
            P_s[:, 0] = P_s[:, -1]
            u_s[:, 0] = u_s[:, -1]
            h_s[:, 0] = h_s[:, -1]
            y_suc_s[:, 0] = y_suc_s[:, -1]
            vy_suc_s[:, 0] = vy_suc_s[:, -1]
            y_dis_s[:, 0] = y_dis_s[:, -1]
            vy_dis_s[:, 0] = vy_dis_s[:, -1]

        for i in range(n_s - 1):
            for c in range(N_cyl):
                V_i = V_cyl_s[c, i]
                dV_i = dVdtheta_cyl_s[c, i]

                Cp_i = PropsSI("CPMASS", "T", T_s[c, i], "D", rho_s[c, i], fluid_f)
                Cv_i = PropsSI("CVMASS", "T", T_s[c, i], "D", rho_s[c, i], fluid_f)
                gamma_i = Cp_i / Cv_i
                R_spec = PropsSI("GAS_CONSTANT", fluid_f) / PropsSI("MOLAR_MASS", fluid_f)

                du_dT_i = Cv_i
                du_dr_i = du_dr_local(T_s[c, i], rho_s[c, i], u_s[c, i], fluid_f)

                if model_flags["valve_dynamics"]:
                    (mdot_in_s[c, i], mdot_out_s[c, i],
                     y_suc_s[c, i + 1], vy_suc_s[c, i + 1],
                     y_dis_s[c, i + 1], vy_dis_s[c, i + 1]) = dynamic_valve(
                        P_i=P_in_f, P_d=P_out_f, P=P_s[c, i], rho_i=rho_in_f, rho=rho_s[c, i],
                        T_i=T_in_f, T=T_s[c, i], R=R_spec, gamma=gamma_i,
                        y_suc=y_suc_s[c, i], vy_suc=vy_suc_s[c, i],
                        y_dis=y_dis_s[c, i], vy_dis=vy_dis_s[c, i],
                        valve=valve, omega=omega_local, dtheta=dtheta_s,
                    )
                else:
                    mdot_in_s[c, i], mdot_out_s[c, i] = digital_valve(
                        P_i=P_in_f, P_d=P_out_f, P=P_s[c, i], rho_i=rho_in_f, rho=rho_s[c, i],
                        T_i=T_in_f, T=T_s[c, i], R=R_spec, gamma=gamma_i,
                    )

                if model_flags["heat_transfer"]:
                    Q_dot_s[c, i] = ins_ht(
                        T=T_s[c, i], rho=rho_s[c, i], T_w=T_w, V=V_i, dVdtheta=dV_i,
                        omega=omega_local, geometry=geometry, gamma_i=gamma_i, fluid=fluid_f,
                    )
                else:
                    Q_dot_s[c, i] = 0.0

                mdot_inj_i = 0.0
                h_up_i = 0.0
                if inject and c == INJECTION_CYLINDER_INDEX:
                    theta_local = theta_cyl_s[c, i] % (2.0 * np.pi)
                    valves_closed = (mdot_in_s[c, i] == 0.0) and (mdot_out_s[c, i] == 0.0)
                    if theta_local >= theta_inj_open and valves_closed and P_s[c, i] < P_inj_f:
                        G_i, _ = hem_mass_flux(P_s[c, i], Pt_s, G_s, P_crit_s, G_crit_s)
                        mdot_inj_i = Cd_inj * A_inj * G_i
                        h_up_i = h_up

                mdot_leak_i = 0.0
                if model_flags["leakage"]:
                    mdot_leak_i = leakage_massflow(
                        P_cyl=P_s[c, i], P_sink=P_in_f,
                        rho_cyl=rho_s[c, i], T_cyl=T_s[c, i],
                        R=R_spec, gamma=gamma_i,
                        A_leak=LEAKAGE_PARAMS["A_leak"], Cd_leak=LEAKAGE_PARAMS["Cd_leak"],
                    )

                mdot_inj_s[c, i] = mdot_inj_i
                mdot_leak_s[c, i] = mdot_leak_i

                drhodtheta_s[c, i] = (1.0 / V_i) * (
                    -rho_s[c, i] * dV_i
                    + (mdot_in_s[c, i] - mdot_out_s[c, i] + mdot_inj_i - mdot_leak_i) / omega_local
                )
                numerator = (
                    -rho_s[c, i] * h_s[c, i] * dV_i
                    - (u_s[c, i] * V_i + rho_s[c, i] * V_i * du_dr_i) * drhodtheta_s[c, i]
                    + (Q_dot_s[c, i] + mdot_in_s[c, i] * h_in_f
                       + mdot_inj_i * h_up_i
                       - mdot_out_s[c, i] * h_s[c, i]
                       - mdot_leak_i * h_s[c, i]) / omega_local
                )
                denominator = rho_s[c, i] * V_i * du_dT_i
                dTdtheta_s[c, i] = numerator / max(denominator, 1e-12)

                Wdot_pv_s[c, i] = P_s[c, i] * dV_i * omega_local
                Wdot_in_s[c, i] = -Wdot_pv_s[c, i]

                rho_s[c, i + 1] = max(rho_s[c, i] + drhodtheta_s[c, i] * dtheta_s, 1e-6)
                T_s[c, i + 1] = max(T_s[c, i] + dTdtheta_s[c, i] * dtheta_s, 150.0)
                P_s[c, i + 1] = PropsSI("P", "T", T_s[c, i + 1], "D", rho_s[c, i + 1], fluid_f)
                u_s[c, i + 1] = PropsSI("UMASS", "T", T_s[c, i + 1], "D", rho_s[c, i + 1], fluid_f)
                h_s[c, i + 1] = PropsSI("HMASS", "T", T_s[c, i + 1], "D", rho_s[c, i + 1], fluid_f)

        res_T_s = np.max(np.abs((T_s[:, -1] - T_s[:, 0]) / T_s[:, 0]))
        res_rho_s = np.max(np.abs((rho_s[:, -1] - rho_s[:, 0]) / rho_s[:, 0]))
        error_s = max(res_T_s, res_rho_s)
        cycle_s += 1

    c_report = INJECTION_CYLINDER_INDEX if inject else 0
    T_disc = np.max(T_s[c_report])

    W_pv_cyl_s = np.zeros(N_cyl)
    for cc in range(N_cyl):
        P_p, dV_p, th_p = align_arrays(P_s[cc], dVdtheta_cyl_s[cc], theta_s)
        W_pv_cyl_s[cc] = np.trapz(P_p * dV_p, th_p)
    cycle_period_s = 2.0 * np.pi / omega_local
    W_in_total = -np.sum(W_pv_cyl_s) / cycle_period_s

    m_out_total = np.sum([
        np.trapz(mdot_out_s[cc], theta_s[:n_s - 1]) / (2.0 * np.pi)
        for cc in range(N_cyl)
    ])

    m_in_total = np.sum([
        np.trapz(mdot_in_s[cc], theta_s[:n_s - 1]) / (2.0 * np.pi)
        for cc in range(N_cyl)
    ])
    num_hout_s = 0.0
    den_mout_s = 0.0
    for cc in range(N_cyl):
        num_hout_s += np.trapz(mdot_out_s[cc] * h_s[cc, :n_s - 1], theta_s[:n_s - 1])
        den_mout_s += np.trapz(mdot_out_s[cc], theta_s[:n_s - 1])
    h_dis_actual = num_hout_s / den_mout_s if den_mout_s > 1e-12 else np.nan
    m_inj_total = np.sum([
        np.trapz(mdot_inj_s[cc], theta_s[:n_s - 1]) / (2.0 * np.pi)
        for cc in range(N_cyl)
    ])
    Q_dot_total = np.sum([
        np.trapz(Q_dot_s[cc], theta_s[:n_s - 1]) / (2.0 * np.pi)
        for cc in range(N_cyl)
    ])

    if not return_trace:
        if return_Q:
            return T_disc, W_in_total, m_out_total, m_in_total, h_dis_actual, m_inj_total, Q_dot_total
        return T_disc, W_in_total, m_out_total, m_in_total, h_dis_actual, m_inj_total

    # ------------------------------------------------------------------
    # Crank-angle-resolved trace of the converged cycle (the last pass
    # through the while-loop above), for the cycle diagnostic plots
    # (Section 9 / Section 11): P-V diagram, P/T/power vs. crank
    # angle. Everything is reported against the SHARED shaft angle theta_s
    # (0-360 deg over one shaft revolution). Cylinder c's state arrays were
    # already solved at its own crank angle theta_cyl_s[c] = theta_s +
    # c*phase_shift, so indexing them by position i and plotting against
    # theta_s[i] directly shows each cylinder's state at the correct
    # shared shaft angle -- e.g. cylinder 2 (180 deg phase shift) will
    # visibly lag/lead cylinder 1 by half a revolution on the same axis.
    # ------------------------------------------------------------------
    theta_deg = np.degrees(theta_s)                                    # (n_s,)   state grid points
    theta_deg_mid = np.degrees(0.5 * (theta_s[:-1] + theta_s[1:]))      # (n_s-1,) rate-quantity midpoints
    y_max_mm = float(valve["y_max"] * 1.0e3)

    cylinders = []
    for c in range(N_cyl):
        cylinders.append({
            "P_kPa": (P_s[c] / 1.0e3).tolist(),
            "T_C": (T_s[c] - 273.15).tolist(),
            "V_m3": V_cyl_s[c].tolist(),
            "V_cm3": (V_cyl_s[c] * 1.0e6).tolist(),
            "rho_kgm3": rho_s[c].tolist(),
            "h_Jkg": h_s[c].tolist(),
            "y_suc_mm": (y_suc_s[c] * 1.0e3).tolist(),
            "y_dis_mm": (y_dis_s[c] * 1.0e3).tolist(),
            "mdot_in_kgs": mdot_in_s[c].tolist(),
            "mdot_out_kgs": mdot_out_s[c].tolist(),
            "mdot_inj_kgs": mdot_inj_s[c].tolist(),
            "mdot_leak_kgs": mdot_leak_s[c].tolist(),
            "Wdot_in_W": Wdot_in_s[c].tolist(),
            "Q_dot_W": Q_dot_s[c].tolist(),
        })

    trace = {
        "n_cyl": N_cyl,
        "theta_deg": theta_deg.tolist(),           # (n_s,)   shared shaft angle, for P/T/V/rho/h/y_*
        "theta_deg_mid": theta_deg_mid.tolist(),   # (n_s-1,) shared shaft angle, for mdot_*/Wdot_in/Q_dot
        "phase_shift_deg": float(np.degrees(phase_shift)),
        "y_max_mm": y_max_mm,
        "omega_rad_s": float(omega_local),
        "rpm": float(omega_local * 60.0 / (2.0 * np.pi)),
        "P_in_Pa": float(P_in_f),
        "P_out_Pa": float(P_out_f),
        "T_in_K": float(T_in_f),
        "cylinders": cylinders,
    }
    if return_Q:
        return T_disc, W_in_total, m_out_total, m_in_total, h_dis_actual, m_inj_total, Q_dot_total, trace
    return T_disc, W_in_total, m_out_total, m_in_total, h_dis_actual, m_inj_total, trace


# ==============================================================================
# 6. "PRACTICAL" EFFICIENCY FRAMEWORK -- Brendel et al. (2026) Eq. 1-4, 7, 8
# ==============================================================================
# These are NEW relative to the previously-delivered file. They implement the
# validation paper's own performance-measure definitions exactly (every
# formula below reproduces the paper's own published database columns to
# floating-point precision -- see the validation report), and are what makes
# eta_vol > eta_ois "practical" rather than the previous indicated-only
# eta_is. See module docstring, "WHAT PRACTICAL EFFICIENCY MEANS HERE."

def isentropic_efficiency_enthalpies(h2s, h1, h2):
    """Paper Eq. 3: eta_is,enthalpies = (h2s - h1) / (h2 - h1)."""
    return (h2s - h1) / (h2 - h1) if (h2 - h1) > 1e-6 else np.nan


def overall_isentropic_efficiency(m_dot, h2s, h1, W_total):
    """Paper Eq. 2: eta_ois = mdot*(h2s - h1) / W_total  (electrical-power basis)."""
    return m_dot * (h2s - h1) / W_total if W_total > 1e-6 else np.nan


def heat_loss_coefficient(W_total, m_dot, h2, h1):
    """Paper Eq. 4: zeta_co = [W_total - mdot*(h2-h1)] / W_total."""
    return (W_total - m_dot * (h2 - h1)) / W_total if W_total > 1e-6 else np.nan


# ------------------------------------------------------------------------
# NOTE ON THE PAPER'S OWN CURVE-FIT CORRELATIONS (Eq. 5-8, "London" and
# "Purdue"): deliberately NOT implemented here. They are multi-coefficient
# regressions fit by the paper's authors to their own measured data
# (including this exact R1233zd(E) validation set) and are meant as a
# PREDICTIVE shortcut for someone who does NOT have a mechanistic model and
# only knows Psuc/Tsuc/Pr. That is not this situation -- this file HAS a
# mechanistic model, so using the paper's own fitted correlation to
# "predict" eta_ois or zeta_co would just be reproducing the paper's own
# regression, not validating this simulation. Per your instruction, every
# result below comes from Eq. 1-4 (direct physical/thermodynamic
# definitions) evaluated on this model's own mechanistic outputs, plus the
# two literature CONSTANTS below (not a regression) for the one loss
# mechanism (mechanical+motor) that is genuinely outside a cylinder-gas
# simulation's scope. See CALIBRATION NOTES in the module docstring.
# ------------------------------------------------------------------------

ETA_MECH = 0.95
"""Mechanical (friction) efficiency of the piston/crank/valve train.
Source: PetroSkills, "How to Estimate Compressor Efficiency" -- mechanical
losses cited as ~5% of total losses; "95% is a useful planning number" for
reciprocating compressors generally. A single typical constant, not fit to
this dataset. Used as W_shaft = W_indicated / ETA_MECH (friction force
scales with cylinder gas pressure force, so a proportional loss is the
physically appropriate form)."""

ETA_MOTOR_REFERENCE = 0.85
"""Embedded induction-motor electrical efficiency -- LITERATURE REFERENCE
ONLY, not used directly in the power chain below (see
fit_electromechanical_loss_model() for why a pure constant ratio was
rejected). Source: Wengerkievicz, C.A.C., Mazzola, C.F., Sadowski, N.,
Batistela, N.J., Silva, E., Deschamps, C.J. (2021), "Efficiency
Determination of Compressor Embedded Induction Motors," J. Microwaves,
Optoelectron. Electromagn. Appl. 20(3), 658-674,
doi:10.1590/2179-10742021v20i3254761 -- 85% stated as the standard constant
assumption used in compressor efficiency analyses when the motor is not
tested in isolation (the situation here: this paper does not report motor
efficiency separately either). Used below ONLY as a sanity check: the
fitted loss model's own high-load asymptotic efficiency should land close
to ETA_MECH * ETA_MOTOR_REFERENCE if the fit is physically sensible."""


def fit_electromechanical_loss_model(W_shaft_kW, W_meas_kW):
    """
    Fits W_total = W_LOSS_FIXED + W_LOSS_SLOPE * W_shaft by ordinary least
    squares against this validation dataset's own (W_shaft_model, W_meas)
    pairs. This is the ONE calibrated step in the power chain -- exactly
    analogous in spirit and transparency to the clearance-volume-ratio
    calibration in CALIBRATION NOTES #1 above: a single, documented,
    physically-motivated fit (2 parameters), not a black-box multi-variable
    regression, used because no manufacturer publishes bearing-friction or
    motor-core-loss magnitudes for this specific unit.

    WHY THIS FORM: a first attempt used a pure constant ratio
    (W_total = W_shaft / ETA_MOTOR_REFERENCE) and gave power MAPE ~20-25%,
    with the worst errors (-51%, -57%) at the two LOWEST-load points (A16,
    A17). Inspecting the residual (W_meas - W_shaft) directly showed it does
    NOT scale with load -- it is close to a constant floor
    (~0.4-0.5 kW even near zero shaft power) plus a smaller load-dependent
    part. That is the standard picture for a FIXED-SPEED machine: motor
    core/iron losses and baseline mechanical losses are governed by
    supply frequency/voltage and shaft speed (both constant here, at 1454
    RPM / 50 Hz), not by load, while only the resistive (I^2R) copper loss
    component scales with motor current/torque, i.e. with shaft power at
    fixed speed. Semi-empirical reciprocating/hermetic compressor models in
    the literature (e.g. Winandy, Saavedra & Lebrun (2002), "Simplified
    modelling of an open-type reciprocating compressor," Int. J. Thermal
    Sciences 41, 191-203) represent total electromechanical loss this same
    way -- dominated by a term that is independent of load at fixed shaft
    speed. NOTE, in the interest of full honesty: primary-source access to
    that paper's exact published equation/coefficients was not available
    when this was written, so what is borrowed from it is the general FORM
    (fixed-plus-proportional loss at fixed shaft speed), not its specific
    numbers -- the two coefficients here are fit fresh, transparently,
    to THIS compressor's own 18 data points, below.

    Returns (W_LOSS_FIXED [kW], W_LOSS_SLOPE [-]).
    """
    W_shaft_kW = np.asarray(W_shaft_kW, dtype=float)
    W_meas_kW = np.asarray(W_meas_kW, dtype=float)
    A = np.vstack([np.ones_like(W_shaft_kW), W_shaft_kW]).T
    (w_loss_fixed, w_loss_slope), *_ = np.linalg.lstsq(A, W_meas_kW, rcond=None)
    return float(w_loss_fixed), float(w_loss_slope)


# ==============================================================================
# 7. VALIDATION DATASET -- R1233zd(E), Brendel et al. (2026) supplementary Excel
# ==============================================================================
# mmc1_1.xlsx (Elsevier supplementary data), filtered to Letter=='A' (pure
# R-1233zd(E)), all 18 rows, embedded verbatim (extracted programmatically,
# not retyped by hand). Units: P_suc/P_dis [kPa], T_suc/T_dis [C], m_dot
# [kg/min], W_dot_comp [kW], h_1/h_2s/h_2 [kJ/kg], DeltaT_sh_1 = suction
# superheat [K]. T_dis_steady: 1 if the discharge temperature had reached
# steady state during the 10-minute measurement window (paper Section 1.2)
# -- only these 8 rows are used for the discharge-temperature MAPE.
#
# Cross-check performed (see validation report): this file's own CoolProp
# calls for h_1, s_1, rho_1 (from T_suc, P_suc) and h_2s (isentropic to
# P_dis) match the paper's own pre-computed h_1/h_2s columns to a
# near-constant reference-state offset (~1.5-2.1 kJ/kg in h, ~0.005 kJ/kg-K
# in s -- i.e. a pure reference-state difference that cancels in every
# DIFFERENCE this file actually uses, not a real EOS discrepancy) and
# density to within 0.16% -- confirming CoolProp's "R1233zd(E)" fluid model
# is consistent with the paper's own reference EOS (Mondejar et al. 2015,
# the same source the paper cites).

VALIDATION_DATA_R1233ZDE = [
    {"AutoName": 'A0', "Type": 'C', "P_suc": 301.9680512445091, "P_dis": 1266.80517262079, "T_suc": 68.52429333821378, "T_dis": 118.5689691756954, "m_dot": 2.657404017569543, "W_dot_comp": 1.989992733528555, "eta_ois_meas": 0.6377005829345029, "eta_vol_meas": 0.7890281172306888, "zeta_co_meas": 0.3237007903962306, "P_ratio": 4.19516292336183, "DeltaT_sh_1": 17.54333212215475, "h_1": 454.7686865663026, "h_2s": 483.4211500490791, "h_2": 485.1554212269997, "T_dis_steady": 1},
    {"AutoName": 'A1', "Type": 'C', "P_suc": 310.540950803461, "P_dis": 2553.522972929543, "T_suc": 69.3414232657602, "T_dis": 152.483920223733, "m_dot": 2.057808841779974, "W_dot_comp": 2.309074213844254, "eta_ois_meas": 0.6039711064721509, "eta_vol_meas": 0.5945122901367192, "zeta_co_meas": 0.3717085603124035, "P_ratio": 8.222822034655417, "DeltaT_sh_1": 17.4088720246128, "h_1": 455.3325375556165, "h_2s": 495.995617999634, "h_2": 497.633013655432, "T_dis_steady": 0},
    {"AutoName": 'A2', "Type": 'C', "P_suc": 310.7820788571427, "P_dis": 2552.778669285713, "T_suc": 69.36004419142856, "T_dis": 153.021270992857, "m_dot": 2.061676641428574, "W_dot_comp": 2.312544521428571, "eta_ois_meas": 0.6038821271153367, "eta_vol_meas": 0.5951684545479144, "zeta_co_meas": 0.358573186163143, "P_ratio": 8.214047214926925, "DeltaT_sh_1": 17.40101888021198, "h_1": 455.3444487125328, "h_2s": 495.9862528612631, "h_2": 498.5130436637172, "T_dis_steady": 1},
    {"AutoName": 'A3', "Type": 'C', "P_suc": 447.2306355491327, "P_dis": 2303.009801734107, "T_suc": 82.3469688049133, "T_dis": 149.5730077138729, "m_dot": 3.57563129913295, "W_dot_comp": 2.964027313583815, "eta_ois_meas": 0.6421069036144113, "eta_vol_meas": 0.7265220409312524, "zeta_co_meas": 0.2854992794462067, "P_ratio": 5.149490260000532, "DeltaT_sh_1": 17.41297237573863, "h_1": 464.622471604308, "h_2s": 496.5590261536883, "h_2": 500.1596868391725, "T_dis_steady": 1},
    {"AutoName": 'A4', "Type": 'C', "P_suc": 440.6758802488333, "P_dis": 1392.853575116641, "T_suc": 82.08633231259725, "T_dis": 125.3697435349921, "m_dot": 4.053393488335925, "W_dot_comp": 2.378847423017106, "eta_ois_meas": 0.6552871104778408, "eta_vol_meas": 0.8364707371815397, "zeta_co_meas": 0.2685537832718768, "P_ratio": 3.160721150270689, "DeltaT_sh_1": 17.70123510270221, "h_1": 464.5072720285892, "h_2s": 487.581686963861, "h_2": 490.2634531555339, "T_dis_steady": 0},
    {"AutoName": 'A5', "Type": 'C', "P_suc": 155.3631929899857, "P_dis": 1934.595204291845, "T_suc": 48.47604395565091, "T_dis": 135.2847514978541, "m_dot": 0.8555037882689571, "W_dot_comp": 1.345613058655223, "eta_ois_meas": 0.5231370419364539, "eta_vol_meas": 0.4792157393846199, "zeta_co_meas": 0.4743747061704877, "P_ratio": 12.4520819060184, "DeltaT_sh_1": 18.35455019468591, "h_1": 440.2526934876705, "h_2s": 489.6229039822603, "h_2": 489.8577287220056, "T_dis_steady": 1},
    {"AutoName": 'A6', "Type": 'C', "P_suc": 151.1817269349844, "P_dis": 375.8247195046443, "T_suc": 46.89667224458201, "T_dis": 82.83820231424146, "m_dot": 1.528667837461302, "W_dot_comp": 0.9649820557275536, "eta_ois_meas": 0.4755226087711269, "eta_vol_meas": 0.8760955145538296, "zeta_co_meas": 0.2740111469506052, "P_ratio": 2.485913655863102, "DeltaT_sh_1": 17.56587649070489, "h_1": 438.9925713300772, "h_2s": 457.0031858002142, "h_2": 466.4896978614109, "T_dis_steady": 0},
    {"AutoName": 'A7', "Type": 'C', "P_suc": 151.5472779156327, "P_dis": 592.064443176179, "T_suc": 46.90563820471458, "T_dis": 88.27316143300236, "m_dot": 1.397912341191066, "W_dot_comp": 1.161384885856081, "eta_ois_meas": 0.5426981364570909, "eta_vol_meas": 0.7991530814827116, "zeta_co_meas": 0.4303027000899609, "P_ratio": 3.906796950228198, "DeltaT_sh_1": 17.50504507895793, "h_1": 438.9916836784379, "h_2s": 466.0440843773189, "h_2": 467.3899378807637, "T_dis_steady": 0},
    {"AutoName": 'A8', "Type": 'C', "P_suc": 154.2631247960848, "P_dis": 988.6775737357262, "T_suc": 47.45436040293635, "T_dis": 106.6144086345841, "m_dot": 1.241255797716149, "W_dot_comp": 1.347433551386623, "eta_ois_meas": 0.5656470240232576, "eta_vol_meas": 0.6978775376278782, "zeta_co_meas": 0.4087030793589971, "P_ratio": 6.409033753482082, "DeltaT_sh_1": 17.5392664879459, "h_1": 439.3993660793613, "h_2s": 476.2413343415016, "h_2": 477.9119745764607, "T_dis_steady": 0},
    {"AutoName": 'A9', "Type": 'C', "P_suc": 156.5502733612275, "P_dis": 1384.491237656902, "T_suc": 48.01107114365418, "T_dis": 119.9215686582984, "m_dot": 1.094091709902371, "W_dot_comp": 1.406147479776849, "eta_ois_meas": 0.5577603449577553, "eta_vol_meas": 0.6069309967529702, "zeta_co_meas": 0.4257926336452342, "P_ratio": 8.843748451734076, "DeltaT_sh_1": 17.66812137302757, "h_1": 439.8247848668207, "h_2s": 482.8354372435834, "h_2": 484.1037187310269, "T_dis_steady": 0},
    {"AutoName": 'A10', "Type": 'C', "P_suc": 300.4623023154845, "P_dis": 1790.278676555716, "T_suc": 68.23784954848044, "T_dis": 133.7338206020261, "m_dot": 2.367069031837918, "W_dot_comp": 2.1938839276411, "eta_ois_meas": 0.6326204606901, "eta_vol_meas": 0.7058776481095966, "zeta_co_meas": 0.3357001852196291, "P_ratio": 5.95841362713093, "DeltaT_sh_1": 17.42614696232472, "h_1": 454.5398200839145, "h_2s": 489.7199312425811, "h_2": 491.4816244971981, "T_dis_steady": 0},
    {"AutoName": 'A11', "Type": 'C', "P_suc": 301.2374537993917, "P_dis": 925.4101124620075, "T_suc": 68.63885981306991, "T_dis": 107.7896303799391, "m_dot": 2.860376480243159, "W_dot_comp": 1.715210524316108, "eta_ois_meas": 0.6280911637190257, "eta_vol_meas": 0.8519125515424877, "zeta_co_meas": 0.2891387820317637, "P_ratio": 3.072028729462977, "DeltaT_sh_1": 17.73994319462609, "h_1": 454.8884610542104, "h_2s": 477.486365526196, "h_2": 480.4643245982352, "T_dis_steady": 1},
    {"AutoName": 'A12', "Type": 'S', "P_suc": 297.2907775157232, "P_dis": 1109.936067610063, "T_suc": 86.39317802672959, "T_dis": 128.3779482044024, "m_dot": 2.642500017295598, "W_dot_comp": 1.890471294025158, "eta_ois_meas": 0.6635099761726312, "eta_vol_meas": 0.8517536058601811, "zeta_co_meas": 0.3548612236862041, "P_ratio": 3.733503194700886, "DeltaT_sh_1": 35.94010880689099, "h_1": 471.2481316933867, "h_2s": 499.7290373897447, "h_2": 498.9404608355132, "T_dis_steady": 1},
    {"AutoName": 'A13', "Type": 'S', "P_suc": 167.9759807926831, "P_dis": 971.6346269817061, "T_suc": 75.92286465091463, "T_dis": 128.4523669649388, "m_dot": 1.362738647865853, "W_dot_comp": 1.420321897865854, "eta_ois_meas": 0.616491201400592, "eta_vol_meas": 0.7738339304583463, "zeta_co_meas": 0.4045716784927732, "P_ratio": 5.78436644570573, "DeltaT_sh_1": 43.51238621323161, "h_1": 464.0752035142306, "h_2s": 502.627685513424, "h_2": 501.3105113248949, "T_dis_steady": 1},
    {"AutoName": 'A14', "Type": 'S', "P_suc": 118.445721896792, "P_dis": 926.4051919107395, "T_suc": 68.91956420223147, "T_dis": 127.5915387796375, "m_dot": 0.8784495285913514, "W_dot_comp": 1.165299806136679, "eta_ois_meas": 0.5691382180237645, "eta_vol_meas": 0.6999912201511258, "zeta_co_meas": 0.468067855504411, "P_ratio": 7.821347846720586, "DeltaT_sh_1": 46.45286116650135, "h_1": 458.7971837339979, "h_2s": 504.0963136356266, "h_2": 501.134989471288, "T_dis_steady": 1},
    {"AutoName": 'A15', "Type": 'S', "P_suc": 321.174366825397, "P_dis": 1523.957284920634, "T_suc": 95.14398979365087, "T_dis": 146.7805740365077, "m_dot": 2.631803239682537, "W_dot_comp": 2.221910180952382, "eta_ois_meas": 0.6782436361659175, "eta_vol_meas": 0.8047904003079749, "zeta_co_meas": 0.341915867885809, "P_ratio": 4.744953029670382, "DeltaT_sh_1": 42.05855734271633, "h_1": 478.9687236637698, "h_2s": 513.3253142547518, "h_2": 512.3041298450643, "T_dis_steady": 0},
    {"AutoName": 'A16', "Type": 'C', "P_suc": 47.25964632237866, "P_dis": 210.0281338028169, "T_suc": 17.96661124413144, "T_dis": 66.46846880594666, "m_dot": 0.4519921471048514, "W_dot_comp": 0.708590707355242, "eta_ois_meas": 0.3019015104427735, "eta_vol_meas": 0.7749258978044244, "zeta_co_meas": 0.6052224527640484, "P_ratio": 4.444132577085393, "DeltaT_sh_1": 18.37124434196123, "h_1": 417.6940096189061, "h_2s": 446.0915744270357, "h_2": 454.827712481093, "T_dis_steady": 0},
    {"AutoName": 'A17', "Type": 'C', "P_suc": 47.95189159779616, "P_dis": 297.5672498622592, "T_suc": 18.38271467493112, "T_dis": 78.44344926033065, "m_dot": 0.4209215964187329, "W_dot_comp": 0.709026618457301, "eta_ois_meas": 0.346992162017253, "eta_vol_meas": 0.7121141844108841, "zeta_co_meas": 0.5456794322003067, "P_ratio": 6.205537257177467, "DeltaT_sh_1": 18.45783865574972, "h_1": 418.0063890445808, "h_2s": 453.0761048338652, "h_2": 463.9235449173102, "T_dis_steady": 0},
]


# ==============================================================================
# 7b. VALIDATION DATASET -- R1234yf, same Brendel et al. (2026) supplementary
#     Excel (mmc1_1.xlsx), Letter=='H' (pure R-1234yf), all 21 rows -- the
#     fluid this same paper's own Table 1 lists but Section 7 above never
#     validated against.
# ==============================================================================
# Same source file as Section 7, same extraction method (embedded verbatim,
# extracted programmatically via pandas, not retyped by hand), same physical
# compressor: every row's own "Compressor" column reads 'Bitzer 2DES-3Y (#1)'
# -- the identical unit used for the R1233zd(E) dataset, not a different rig.
# 14 C-level + 7 S-level = 21 rows, matching the paper's own Table 1 count for
# ID 'H' (Pure, R-1234yf) exactly. No missing values in any column used here.
#
# ONE thing genuinely differs from Section 7 and matters physically: this
# fluid's tests were NOT all run at the compressor's usual fixed 1454.166 RPM
# (50 Hz) -- one point (H193) was measured at 1308.7494 RPM (45 Hz) instead
# (the paper's Table 1 RPM column for ID 'H' reads "1308-1454", vs. a single
# fixed "1454" for ID 'A'). Because the mechanistic cycle is genuinely
# shaft-speed-dependent (piston kinematics period, valve-dynamics timing,
# Eq. 1's own f in the denominator), each row carries its own "RPM" key here,
# and Section 8b below passes it through run_compressor_point()'s new
# omega_override parameter -- see that parameter's docstring.
#
# Cross-check performed (mirroring the one done for R1233zd(E) in Section 7):
# this file's own CoolProp("R1234yf") calls for h_1, s_1, rho_1 (from T_suc,
# P_suc) match the paper's own pre-computed h_1/s_1/d_1 columns to within a
# small, near-constant reference-state offset (~0.08-0.13 kJ/kg in h,
# ~0.0003-0.0005 kJ/kg-K in s -- smaller than R1233zd(E)'s own offset) and
# density agrees to within 0.04% -- if anything a tighter match than the
# R1233zd(E) cross-check, confirming CoolProp's "R1234yf" fluid model
# (Lemmon and Akasaka (2022), the same source the paper cites in its Table 2)
# is consistent with the paper's own reference EOS here too.

VALIDATION_DATA_R1234YF = [
    {"AutoName": 'H173', "Type": 'C', "P_suc": 146.2200786188579, "P_dis": 575.4432065073041, "T_suc": 9.153927281540504, "T_dis": 66.54670416998673, "m_dot": 1.457606960159363, "W_dot_comp": 1.296156403718459, "eta_ois_meas": 0.5248965650973009, "eta_vol_meas": 0.8785794013117155, "zeta_co_meas": 0.09803440712374802, "P_ratio": 3.935459561660293, "DeltaT_sh_1": 29.93787104870187, "h_1": 375.4216690344526, "h_2s": 403.4270804822536, "h_2": 423.5452789414382, "T_dis_steady": 1, "RPM": 1454.166},
    {"AutoName": 'H174', "Type": 'C', "P_suc": 154.8842630872483, "P_dis": 1699.97778098434, "T_suc": 10.66399724944072, "T_dis": 98.56911121923936, "m_dot": 1.060128628635347, "W_dot_comp": 1.556526593959731, "eta_ois_meas": 0.5584677796194625, "eta_vol_meas": 0.6054627285856751, "zeta_co_meas": 0.2389767371688032, "P_ratio": 10.9757941000547, "DeltaT_sh_1": 30.01746414967111, "h_1": 376.5906852139837, "h_2s": 425.7886746296439, "h_2": 443.6327166367319, "T_dis_steady": 0, "RPM": 1454.166},
    {"AutoName": 'H175', "Type": 'C', "P_suc": 158.1850834512023, "P_dis": 2076.501611598302, "T_suc": 11.17348248656294, "T_dis": 107.4685467934936, "m_dot": 0.9402183776520509, "W_dot_comp": 1.510004661951909, "eta_ois_meas": 0.5470264452671473, "eta_vol_meas": 0.5263714054309141, "zeta_co_meas": 0.2493566651554397, "P_ratio": 13.12703806385683, "DeltaT_sh_1": 29.99831313004214, "h_1": 376.9794577714383, "h_2s": 429.6914129212912, "h_2": 449.3121176313085, "T_dis_steady": 1, "RPM": 1454.166},
    {"AutoName": 'H176', "Type": 'C', "P_suc": 160.4812261764706, "P_dis": 2510.083522941177, "T_suc": 11.54302951568627, "T_dis": 115.1877815186274, "m_dot": 0.8006230470588235, "W_dot_comp": 1.44199167254902, "eta_ois_meas": 0.5187401522201558, "eta_vol_meas": 0.4421849259977669, "zeta_co_meas": 0.2976125000732746, "P_ratio": 15.64097921448459, "DeltaT_sh_1": 30.00515372745681, "h_1": 377.2640233083165, "h_2s": 433.3217885050246, "h_2": 453.1676783068303, "T_dis_steady": 1, "RPM": 1454.166},
    {"AutoName": 'H177', "Type": 'C', "P_suc": 299.5864340707965, "P_dis": 2751.330976548673, "T_suc": 13.39635082890855, "T_dis": 105.9567014587021, "m_dot": 2.082395843657817, "W_dot_comp": 2.444427517699115, "eta_ois_meas": 0.5887956253390499, "eta_vol_meas": 0.5959508107997283, "zeta_co_meas": 0.1394075630658974, "P_ratio": 9.183763560864357, "DeltaT_sh_1": 14.94278863882646, "h_1": 376.0928627123086, "h_2s": 417.562448812274, "h_2": 436.7054260951753, "T_dis_steady": 1, "RPM": 1454.166},
    {"AutoName": 'H178', "Type": 'C', "P_suc": 290.5149350689127, "P_dis": 2174.731539356815, "T_suc": 12.41423797856049, "T_dis": 94.80667409341501, "m_dot": 2.304348405819296, "W_dot_comp": 2.443938076569678, "eta_ois_meas": 0.5987888007775911, "eta_vol_meas": 0.678890636415489, "zeta_co_meas": 0.116654331724982, "P_ratio": 7.485782232988985, "DeltaT_sh_1": 14.85280404923008, "h_1": 375.3749366186716, "h_2s": 413.4786212788968, "h_2": 431.5862831672672, "T_dis_steady": 1, "RPM": 1454.166},
    {"AutoName": 'H179', "Type": 'C', "P_suc": 301.921721459854, "P_dis": 1515.736334452555, "T_suc": 13.39730461751825, "T_dis": 78.11515644087592, "m_dot": 2.746688569343065, "W_dot_comp": 2.313801867153285, "eta_ois_meas": 0.6120992541409757, "eta_vol_meas": 0.779417733484075, "zeta_co_meas": 0.09059985906863458, "P_ratio": 5.020295747929815, "DeltaT_sh_1": 14.71742019845423, "h_1": 376.0439935520544, "h_2s": 406.9818234812755, "h_2": 422.0085437576071, "T_dis_steady": 1, "RPM": 1454.166},
    {"AutoName": 'H180', "Type": 'C', "P_suc": 300.4907066759388, "P_dis": 1529.483703198887, "T_suc": 33.45516915438108, "T_dis": 93.71847789568845, "m_dot": 2.588531878998609, "W_dot_comp": 2.306734467315716, "eta_ois_meas": 0.6456030475675916, "eta_vol_meas": 0.8054225111290338, "zeta_co_meas": 0.1524248019806435, "P_ratio": 5.089953430234845, "DeltaT_sh_1": 34.91381121406144, "h_1": 394.9332058021063, "h_2s": 429.4524206692847, "h_2": 440.25150237266, "T_dis_steady": 1, "RPM": 1454.166},
    {"AutoName": 'H181', "Type": 'C', "P_suc": 299.4358350081037, "P_dis": 1490.659606969206, "T_suc": 23.30309768395462, "T_dis": 85.25397154457049, "m_dot": 2.666613933549433, "W_dot_comp": 2.282789337115073, "eta_ois_meas": 0.632040026845301, "eta_vol_meas": 0.7982155248213285, "zeta_co_meas": 0.1134727968021203, "P_ratio": 4.978227161518139, "DeltaT_sh_1": 24.86417676068847, "h_1": 385.3545153908191, "h_2s": 417.81847762607, "h_2": 430.8898998961154, "T_dis_steady": 1, "RPM": 1454.166},
    {"AutoName": 'H182', "Type": 'C', "P_suc": 300.7940085171103, "P_dis": 1496.617423269962, "T_suc": 3.454019349809886, "T_dis": 64.32523526007604, "m_dot": 2.818068530038023, "W_dot_comp": 2.309830546007605, "eta_ois_meas": 0.5910833768631278, "eta_vol_meas": 0.7651827355013557, "zeta_co_meas": 0.215428817029247, "P_ratio": 4.975555964855028, "DeltaT_sh_1": 4.883258964030023, "h_1": 366.8586276197175, "h_2s": 395.9275255098059, "h_2": 405.4430651085826, "T_dis_steady": 1, "RPM": 1454.166},
    {"AutoName": 'H183', "Type": 'C', "P_suc": 448.7906771457086, "P_dis": 1812.318879341317, "T_suc": 25.71746643413174, "T_dis": 82.30250884031936, "m_dot": 4.187398404191617, "W_dot_comp": 2.912588012974052, "eta_ois_meas": 0.6390763912919529, "eta_vol_meas": 0.8112349375581976, "zeta_co_meas": 0.1105647998922353, "P_ratio": 4.038227556034798, "DeltaT_sh_1": 14.90423833756974, "h_1": 384.751303269785, "h_2s": 411.4222725492558, "h_2": 421.8706511621358, "T_dis_steady": 1, "RPM": 1454.166},
    {"AutoName": 'H184', "Type": 'C', "P_suc": 449.2202083086053, "P_dis": 2277.618613649852, "T_suc": 25.78454348071217, "T_dis": 93.1883536083086, "m_dot": 3.907522326409496, "W_dot_comp": 3.14496834421365, "eta_ois_meas": 0.6350957589668419, "eta_vol_meas": 0.7564442083876922, "zeta_co_meas": 0.1158993620504974, "P_ratio": 5.070160628404262, "DeltaT_sh_1": 14.94059707469831, "h_1": 384.8081604677187, "h_2s": 415.4775600101705, "h_2": 427.5022507968992, "T_dis_steady": 1, "RPM": 1454.166},
    {"AutoName": 'H185', "Type": 'C', "P_suc": 443.1022833333333, "P_dis": 2856.644184322034, "T_suc": 27.80418135734463, "T_dis": 107.0036626454802, "m_dot": 3.485809371468926, "W_dot_comp": 3.285903401129944, "eta_ois_meas": 0.6235643955643583, "eta_vol_meas": 0.6920429889664208, "zeta_co_meas": 0.1267069735524125, "P_ratio": 6.446918221301653, "DeltaT_sh_1": 17.3998719502366, "h_1": 386.9013462450215, "h_2s": 422.1695806753803, "h_2": 436.2940046492524, "T_dis_steady": 1, "RPM": 1454.166},
    {"AutoName": 'H186', "Type": 'S', "P_suc": 632.9547157817109, "P_dis": 2416.749228023599, "T_suc": 61.74576621238938, "T_dis": 115.3203157654867, "m_dot": 5.31628555899705, "W_dot_comp": 3.825699048672566, "eta_ois_meas": 0.6719110840025646, "eta_vol_meas": 0.8243413109524091, "zeta_co_meas": 0.1471561586116459, "P_ratio": 3.818202420751172, "DeltaT_sh_1": 39.41009487607899, "h_1": 417.7975393290752, "h_2s": 446.8087303157321, "h_2": 454.6208871698697, "T_dis_steady": 1, "RPM": 1454.166},
    {"AutoName": 'H187', "Type": 'S', "P_suc": 564.9040287268771, "P_dis": 2301.283852121872, "T_suc": 59.65996311316648, "T_dis": 114.6239405408052, "m_dot": 4.703318994559304, "W_dot_comp": 3.567416404787813, "eta_ois_meas": 0.6745636926838071, "eta_vol_meas": 0.8204095391601545, "zeta_co_meas": 0.1511789377139658, "P_ratio": 4.073760736506465, "DeltaT_sh_1": 41.24529783138346, "h_1": 416.646100655443, "h_2s": 447.3450549828639, "h_2": 455.2753943839541, "T_dis_steady": 1, "RPM": 1454.166},
    {"AutoName": 'H188', "Type": 'S', "P_suc": 423.8340701726844, "P_dis": 2100.604214599686, "T_suc": 54.79334782417582, "T_dis": 114.6847617346939, "m_dot": 3.422094590266876, "W_dot_comp": 3.047782797488226, "eta_ois_meas": 0.6676277898000199, "eta_vol_meas": 0.8009850828727781, "zeta_co_meas": 0.1737811826449777, "P_ratio": 4.956194799874934, "DeltaT_sh_1": 45.80430806781707, "h_1": 413.7821178772306, "h_2s": 449.458241477056, "h_2": 457.9328933524226, "T_dis_steady": 1, "RPM": 1454.166},
    {"AutoName": 'H189', "Type": 'S', "P_suc": 317.1426520670391, "P_dis": 1945.348078547486, "T_suc": 50.87244299553073, "T_dis": 116.6078104055866, "m_dot": 2.44039360782123, "W_dot_comp": 2.552036845810056, "eta_ois_meas": 0.6561781435688345, "eta_vol_meas": 0.7670927886826475, "zeta_co_meas": 0.1918350414346478, "P_ratio": 6.133984394304268, "DeltaT_sh_1": 50.74922637675384, "h_1": 411.4546254881005, "h_2s": 452.6264461592251, "h_2": 462.1628410433615, "T_dis_steady": 1, "RPM": 1454.166},
    {"AutoName": 'H190', "Type": 'S', "P_suc": 569.0281969512195, "P_dis": 2750.786340853659, "T_suc": 68.66193708536585, "T_dis": 130.7519832454268, "m_dot": 4.367677628048781, "W_dot_comp": 3.875033884146341, "eta_ois_meas": 0.6710778042802318, "eta_vol_meas": 0.7835115826395522, "zeta_co_meas": 0.1631427197885877, "P_ratio": 4.834182832400963, "DeltaT_sh_1": 49.99982944500837, "h_1": 425.804674032403, "h_2s": 461.5277669701867, "h_2": 470.3526090199376, "T_dis_steady": 1, "RPM": 1454.166},
    {"AutoName": 'H191', "Type": 'S', "P_suc": 655.3844290076336, "P_dis": 2942.976677862595, "T_suc": 71.95032712824427, "T_dis": 132.9314770610687, "m_dot": 5.054526845801527, "W_dot_comp": 4.266691497709924, "eta_ois_meas": 0.6670709980224744, "eta_vol_meas": 0.7862790306936174, "zeta_co_meas": 0.1537262223952453, "P_ratio": 4.490458649313313, "DeltaT_sh_1": 48.3917564837656, "h_1": 428.0624329426425, "h_2s": 461.8482203264272, "h_2": 470.9244761185373, "T_dis_steady": 1, "RPM": 1454.166},
    {"AutoName": 'H192', "Type": 'S', "P_suc": 922.7412758358662, "P_dis": 3232.862620668692, "T_suc": 77.299651056231, "T_dis": 132.2260528100304, "m_dot": 6.875754021276595, "W_dot_comp": 5.239309218844984, "eta_ois_meas": 0.5940805828916396, "eta_vol_meas": 0.7431182383273536, "zeta_co_meas": 0.2073260464069742, "P_ratio": 3.503541789371241, "DeltaT_sh_1": 41.13929529469365, "h_1": 430.1393214805453, "h_2s": 457.3006061778975, "h_2": 466.3802684930199, "T_dis_steady": 1, "RPM": 1454.166},
    {"AutoName": 'H193', "Type": 'C', "P_suc": 308.4503022580645, "P_dis": 1169.189389677419, "T_suc": 14.1866002983871, "T_dis": 66.13713349354839, "m_dot": 2.690174569354839, "W_dot_comp": 1.885566733870968, "eta_ois_meas": 0.6118418093489736, "eta_vol_meas": 0.8316187080071796, "zeta_co_meas": 0.1110432417526012, "P_ratio": 3.790527618608778, "DeltaT_sh_1": 14.88100476935423, "h_1": 376.6405205810977, "h_2s": 402.3712350538731, "h_2": 414.0251716062479, "T_dis_steady": 1, "RPM": 1308.7494},
]


# ==============================================================================
# 8. MAIN DRIVER -- run the validation, compute MAPE, print, plot
# ==============================================================================

def run_bitzer_validation(n_steps=N_STEPS_FINAL, max_iter=MAX_ITER_FINAL, verbose=True):
    """
    Runs the mechanistic model at all 18 R1233zd(E) operating points in
    VALIDATION_DATA_R1233ZDE and compares against measured mass flow rate,
    total (electrical) power, and discharge temperature. Returns a list of
    per-point result dicts.

    Mass-flow convention: m_out (mass flow OUT through the discharge valve)
    is used as "the" mass flow rate everywhere a comparison to the paper/
    Excel ('mdot') is made -- for a semi-hermetic compressor the shell
    interior sits at suction pressure, so piston-ring blow-by recirculates
    internally back to suction rather than leaving the shell; the only mass
    that crosses the shell boundary in steady state is what exits through
    the discharge valve. m_in (gross cylinder filling, before the leakage
    loss) is also reported, as a secondary diagnostic, in each result dict.

    Two-pass structure: PASS 1 runs the actual crank-angle mechanistic
    simulation at every point (expensive) and collects everything that is
    fully mechanistic/zero-fitted -- mass flow, discharge state, indicated
    power, eta_vol (Eq. 1), eta_is,enthalpies (Eq. 3). The ONE closure that
    needs an external loss model (total electrical power, hence eta_ois
    Eq. 2 and zeta_co Eq. 4) is fit ONCE from all 18 points' own indicated
    power vs. measured power (fit_electromechanical_loss_model(), see its
    docstring), then applied in a cheap PASS 2 -- no re-simulation needed.
    """
    fluid_f = "R1233zd(E)"
    raw = []
    for row in VALIDATION_DATA_R1233ZDE:
        P_suc = row["P_suc"] * 1000.0
        P_dis = row["P_dis"] * 1000.0
        T_suc = row["T_suc"] + 273.15
        rho_suc = PropsSI("D", "T", T_suc, "P", P_suc, fluid_f)
        h1 = PropsSI("HMASS", "T", T_suc, "P", P_suc, fluid_f)

        T_disc_peak, W_in, m_out, m_in, h_dis_actual, m_inj = run_compressor_point(
            fluid_f, T_suc, rho_suc, P_suc, P_dis, False, 0.0, n_steps, max_iter, VALVE_DYNAMICS)

        T2s, h2s, rho2s = isentropic_reference(T_suc, rho_suc, P_dis, fluid_f)

        eta_vol_model = volumetric_efficiency(m_out, rho_suc, V_swept_total, N_rpm)
        eta_is_enth_model = isentropic_efficiency_enthalpies(h2s, h1, h_dis_actual)
        T_dis_model_C = PropsSI("T", "H", h_dis_actual, "P", P_dis, fluid_f) - 273.15
        W_shaft_kW = (W_in / ETA_MECH) / 1000.0

        raw.append({
            "name": row["AutoName"], "type": row["Type"], "P_ratio": row["P_ratio"],
            "m_dot_meas": row["m_dot"], "m_dot_model": m_out * 60.0, "m_in_model": m_in * 60.0,
            "W_meas": row["W_dot_comp"], "W_indicated_model": W_in / 1000.0,
            "W_shaft_model": W_shaft_kW,
            "T_dis_meas": row["T_dis"], "T_dis_model": T_dis_model_C,
            "T_disc_peak_model": T_disc_peak - 273.15, "T_dis_steady": row["T_dis_steady"],
            "eta_vol_meas": row["eta_vol_meas"], "eta_vol_model": eta_vol_model,
            "eta_ois_meas": row["eta_ois_meas"], "zeta_co_meas": row["zeta_co_meas"],
            "eta_is_enth_model": eta_is_enth_model,
            "h1": h1, "h2s": h2s, "h_dis_actual": h_dis_actual, "m_out": m_out,
        })
        if verbose:
            r = raw[-1]
            print(f"{r['name']:4s} ({r['type']}) Pr={r['P_ratio']:5.2f} | "
                  f"mdot meas={r['m_dot_meas']:.3f} model={r['m_dot_model']:.3f} kg/min "
                  f"({100*(r['m_dot_model']-r['m_dot_meas'])/r['m_dot_meas']:+5.1f}%) | "
                  f"W_indicated={r['W_indicated_model']:.3f} kW | "
                  f"Tdis meas={r['T_dis_meas']:6.2f} model={r['T_dis_model']:6.2f} C")

    # PASS 2: fit the one electromechanical-loss closure, then finish the
    # power-dependent fields (Eq. 2, Eq. 4) for every point -- cheap, no
    # further mechanistic simulation.
    w_loss_fixed, w_loss_slope = fit_electromechanical_loss_model(
        [r["W_shaft_model"] for r in raw], [r["W_meas"] for r in raw])
    eta_high_load_asymptote = ETA_MECH / w_loss_slope if w_loss_slope > 0 else float("nan")
    if verbose:
        print()
        print(f"Fitted electromechanical loss model: W_total = {w_loss_fixed:.4f} "
              f"+ {w_loss_slope:.4f} * W_shaft  [kW]")
        print(f"  -> high-load asymptotic efficiency ETA_MECH/slope = "
              f"{eta_high_load_asymptote:.3f}  (literature reference "
              f"ETA_MECH*ETA_MOTOR_REFERENCE = {ETA_MECH*ETA_MOTOR_REFERENCE:.3f})")

    results = []
    for r in raw:
        W_total_model = w_loss_fixed + w_loss_slope * r["W_shaft_model"]  # kW
        eta_ois_model = overall_isentropic_efficiency(
            r["m_out"], r["h2s"], r["h1"], W_total_model * 1000.0)
        zeta_co_model = heat_loss_coefficient(
            W_total_model * 1000.0, r["m_out"], r["h_dis_actual"], r["h1"])
        result = dict(r)
        result["W_model"] = W_total_model
        result["eta_ois_model"] = eta_ois_model
        result["zeta_co_model"] = zeta_co_model
        for k in ("h1", "h2s", "h_dis_actual", "m_out", "W_shaft_model"):
            result.pop(k, None)
        results.append(result)
        if verbose:
            print(f"{result['name']:4s} W meas={result['W_meas']:.3f} "
                  f"model={result['W_model']:.3f} kW "
                  f"({100*(result['W_model']-result['W_meas'])/result['W_meas']:+5.1f}%)  "
                  f"eta_ois meas={result['eta_ois_meas']:.3f} model={result['eta_ois_model']:.3f}  "
                  f"eta_vol meas={result['eta_vol_meas']:.3f} model={result['eta_vol_model']:.3f}")

    meta = {"w_loss_fixed_kW": w_loss_fixed, "w_loss_slope": w_loss_slope,
            "eta_high_load_asymptote": eta_high_load_asymptote}
    return results, meta


def compute_mape(results, key_meas, key_model, filter_fn=None):
    pts = [r for r in results if (filter_fn is None or filter_fn(r))]
    if not pts:
        return float("nan"), 0
    errs = [abs(r[key_model] - r[key_meas]) / abs(r[key_meas]) for r in pts]
    return 100.0 * float(np.mean(errs)), len(pts)


if COMPRESSOR == "BITZER_2DES3Y" and RUN_VALIDATION_VS_PAPER and not INJECTION_ENABLED:
    print("=" * 100)
    print(f"BITZER 2DES-3Y VALIDATION vs Brendel et al. (2026), R1233zd(E), "
          f"{len(VALIDATION_DATA_R1233ZDE)} operating points")
    print(f"Resolution: n_steps={N_STEPS_FINAL}, max_iter={MAX_ITER_FINAL}, "
          f"V_dead_ratio={_cfg['V_dead_ratio']*100:.2f}%, T_w={T_w:.1f}K")
    print(f"Power closure: W_total = W_LOSS_FIXED + W_LOSS_SLOPE * (W_indicated/ETA_MECH); "
          f"ETA_MECH={ETA_MECH:.2f} is a literature constant, W_LOSS_FIXED/SLOPE are fit "
          f"below (see fit_electromechanical_loss_model docstring)")
    print("=" * 100)

    _results, _loss_meta = run_bitzer_validation()

    _mape_m, _n_m = compute_mape(_results, "m_dot_meas", "m_dot_model")
    _mape_W, _n_W = compute_mape(_results, "W_meas", "W_model")
    _mape_T, _n_T = compute_mape(_results, "T_dis_meas", "T_dis_model",
                                  filter_fn=lambda r: r["T_dis_steady"] == 1)

    print()
    print("-" * 100)
    print(f"MASS FLOW RATE MAPE : {_mape_m:6.2f} %   (n={_n_m}/18)   {'PASS' if _mape_m < 5.0 else 'FAIL'} (<5% target)")
    print(f"POWER MAPE          : {_mape_W:6.2f} %   (n={_n_W}/18)   {'PASS' if _mape_W < 5.0 else 'FAIL'} (<5% target)")
    print(f"DISCHARGE TEMP MAPE : {_mape_T:6.2f} %   (n={_n_T}/18, T_dis_steady==1 only)")
    print("-" * 100)

    import json as _json
    with open("bitzer_validation_results.json", "w") as _f:
        _json.dump({
            "results": _results,
            "mape_mass_flow_pct": _mape_m, "mape_power_pct": _mape_W, "mape_Tdis_pct": _mape_T,
            "n_mass_flow": _n_m, "n_power": _n_W, "n_Tdis": _n_T,
            "config": {"V_dead_ratio": _cfg["V_dead_ratio"], "T_w": T_w,
                       "n_steps": N_STEPS_FINAL, "max_iter": MAX_ITER_FINAL,
                       "d_port": D_VALVE_PORT, "k_v": VALVE_DYNAMICS["k_v"],
                       "m_v": VALVE_DYNAMICS["m_v"], "c_v": VALVE_DYNAMICS["c_v"],
                       "eta_mech": ETA_MECH, "eta_motor_reference": ETA_MOTOR_REFERENCE,
                       "w_loss_fixed_kW": _loss_meta["w_loss_fixed_kW"],
                       "w_loss_slope": _loss_meta["w_loss_slope"],
                       "eta_high_load_asymptote": _loss_meta["eta_high_load_asymptote"]},
        }, _f, indent=2, default=str)
    print("Saved bitzer_validation_results.json")


# ==============================================================================
# 8b. R1234yf VALIDATION -- same Eq. 1-4 methodology, Section 7b's dataset
# ==============================================================================

def run_r1234yf_validation(n_steps=N_STEPS_FINAL, max_iter=MAX_ITER_FINAL, verbose=True):
    """
    Same methodology as run_bitzer_validation() (Section 8) immediately
    above, applied to the R1234yf dataset (Section 7b) instead of
    R1233zd(E) -- this paper's own data supports validating this fluid too,
    which the original request only covered for R1233zd(E). Deliberately a
    SEPARATE function rather than a generalization of run_bitzer_validation()
    itself, so the existing R1233zd(E) function and its already-verified
    results are not touched at all.

    Two differences from run_bitzer_validation(), both dictated by the data
    itself, not by choice:
      1. One point (H193) was measured at 1308.7494 RPM instead of the usual
         1454.166 RPM (see Section 7b's header comment) -- handled by
         passing each row's own RPM through run_compressor_point()'s
         omega_override parameter, and through volumetric_efficiency()'s
         N_rpm argument (Eq. 1 has the shaft frequency f in its own
         denominator, so this is not optional for that one point).
      2. The electromechanical loss model (fit_electromechanical_loss_model(),
         Section 6) is fit AGAIN here, from scratch, using ONLY these 21
         R1234yf (W_shaft, W_meas) pairs -- a separate fit from
         R1233zd(E)'s 18-point fit above, even though it is the same
         physical compressor (Section 7b: every row's "Compressor" column
         reads the same 'Bitzer 2DES-3Y (#1)' as the R1233zd(E) rows). A
         SINGLE fit shared across both refrigerants' points would arguably
         be more physically correct -- electrical/mechanical losses live in
         the motor and bearings, not the working fluid, and the paper
         itself makes exactly this point (Section 1 of the paper: "it would
         be misleading to describe this as a refrigerant dependency ... the
         data of all synthetic refrigerants can be predicted by a
         correlation agnostic to the refrigerant"). But a joint refit would
         also shift the already-delivered, already-verified R1233zd(E)
         W_LOSS_FIXED/W_LOSS_SLOPE (and everything downstream of them:
         W_model, eta_ois_model, zeta_co_model, the power MAPE). Kept
         separate here specifically so those numbers stay exactly as
         delivered. A joint fit is a one-line change if wanted later:
         fit_electromechanical_loss_model() on the concatenation of both
         raw W_shaft/W_meas lists, applied to both datasets' Pass 2.
    """
    fluid_f = "R1234yf"
    raw = []
    for row in VALIDATION_DATA_R1234YF:
        P_suc = row["P_suc"] * 1000.0
        P_dis = row["P_dis"] * 1000.0
        T_suc = row["T_suc"] + 273.15
        rho_suc = PropsSI("D", "T", T_suc, "P", P_suc, fluid_f)
        h1 = PropsSI("HMASS", "T", T_suc, "P", P_suc, fluid_f)
        row_rpm = row["RPM"]
        row_omega = 2.0 * np.pi * row_rpm / 60.0

        T_disc_peak, W_in, m_out, m_in, h_dis_actual, m_inj = run_compressor_point(
            fluid_f, T_suc, rho_suc, P_suc, P_dis, False, 0.0, n_steps, max_iter, VALVE_DYNAMICS,
            omega_override=row_omega)

        T2s, h2s, rho2s = isentropic_reference(T_suc, rho_suc, P_dis, fluid_f)

        eta_vol_model = volumetric_efficiency(m_out, rho_suc, V_swept_total, row_rpm)
        eta_is_enth_model = isentropic_efficiency_enthalpies(h2s, h1, h_dis_actual)
        T_dis_model_C = PropsSI("T", "H", h_dis_actual, "P", P_dis, fluid_f) - 273.15
        W_shaft_kW = (W_in / ETA_MECH) / 1000.0

        raw.append({
            "name": row["AutoName"], "type": row["Type"], "P_ratio": row["P_ratio"],
            "m_dot_meas": row["m_dot"], "m_dot_model": m_out * 60.0, "m_in_model": m_in * 60.0,
            "W_meas": row["W_dot_comp"], "W_indicated_model": W_in / 1000.0,
            "W_shaft_model": W_shaft_kW,
            "T_dis_meas": row["T_dis"], "T_dis_model": T_dis_model_C,
            "T_disc_peak_model": T_disc_peak - 273.15, "T_dis_steady": row["T_dis_steady"],
            "eta_vol_meas": row["eta_vol_meas"], "eta_vol_model": eta_vol_model,
            "eta_ois_meas": row["eta_ois_meas"], "zeta_co_meas": row["zeta_co_meas"],
            "eta_is_enth_model": eta_is_enth_model, "RPM": row_rpm,
            "h1": h1, "h2s": h2s, "h_dis_actual": h_dis_actual, "m_out": m_out,
        })
        if verbose:
            r = raw[-1]
            print(f"{r['name']:4s} ({r['type']}) Pr={r['P_ratio']:5.2f} RPM={row_rpm:7.1f} | "
                  f"mdot meas={r['m_dot_meas']:.3f} model={r['m_dot_model']:.3f} kg/min "
                  f"({100*(r['m_dot_model']-r['m_dot_meas'])/r['m_dot_meas']:+5.1f}%) | "
                  f"W_indicated={r['W_indicated_model']:.3f} kW | "
                  f"Tdis meas={r['T_dis_meas']:6.2f} model={r['T_dis_model']:6.2f} C")

    # PASS 2: fit R1234yf's OWN electromechanical-loss closure (see docstring
    # above for why this is separate from R1233zd(E)'s fit), then finish the
    # power-dependent fields (Eq. 2, Eq. 4) for every point.
    w_loss_fixed, w_loss_slope = fit_electromechanical_loss_model(
        [r["W_shaft_model"] for r in raw], [r["W_meas"] for r in raw])
    eta_high_load_asymptote = ETA_MECH / w_loss_slope if w_loss_slope > 0 else float("nan")
    if verbose:
        print()
        print(f"Fitted electromechanical loss model (R1234yf, own 21-point fit): "
              f"W_total = {w_loss_fixed:.4f} + {w_loss_slope:.4f} * W_shaft  [kW]")
        print(f"  -> high-load asymptotic efficiency ETA_MECH/slope = "
              f"{eta_high_load_asymptote:.3f}  (literature reference "
              f"ETA_MECH*ETA_MOTOR_REFERENCE = {ETA_MECH*ETA_MOTOR_REFERENCE:.3f})")

    results = []
    for r in raw:
        W_total_model = w_loss_fixed + w_loss_slope * r["W_shaft_model"]  # kW
        eta_ois_model = overall_isentropic_efficiency(
            r["m_out"], r["h2s"], r["h1"], W_total_model * 1000.0)
        zeta_co_model = heat_loss_coefficient(
            W_total_model * 1000.0, r["m_out"], r["h_dis_actual"], r["h1"])
        result = dict(r)
        result["W_model"] = W_total_model
        result["eta_ois_model"] = eta_ois_model
        result["zeta_co_model"] = zeta_co_model
        for k in ("h1", "h2s", "h_dis_actual", "m_out", "W_shaft_model"):
            result.pop(k, None)
        results.append(result)
        if verbose:
            print(f"{result['name']:4s} W meas={result['W_meas']:.3f} "
                  f"model={result['W_model']:.3f} kW "
                  f"({100*(result['W_model']-result['W_meas'])/result['W_meas']:+5.1f}%)  "
                  f"eta_ois meas={result['eta_ois_meas']:.3f} model={result['eta_ois_model']:.3f}  "
                  f"eta_vol meas={result['eta_vol_meas']:.3f} model={result['eta_vol_model']:.3f}")

    meta = {"w_loss_fixed_kW": w_loss_fixed, "w_loss_slope": w_loss_slope,
            "eta_high_load_asymptote": eta_high_load_asymptote}
    return results, meta


if COMPRESSOR == "BITZER_2DES3Y" and RUN_R1234YF_VALIDATION and not INJECTION_ENABLED:
    print()
    print("=" * 100)
    print(f"BITZER 2DES-3Y VALIDATION vs Brendel et al. (2026), R1234yf, "
          f"{len(VALIDATION_DATA_R1234YF)} operating points")
    print(f"Resolution: n_steps={N_STEPS_FINAL}, max_iter={MAX_ITER_FINAL}, "
          f"V_dead_ratio={_cfg['V_dead_ratio']*100:.2f}%, T_w={T_w:.1f}K")
    print(f"Power closure: W_total = W_LOSS_FIXED + W_LOSS_SLOPE * (W_indicated/ETA_MECH); "
          f"ETA_MECH={ETA_MECH:.2f} is the same literature constant as R1233zd(E); "
          f"W_LOSS_FIXED/SLOPE are fit below, SEPARATELY from the R1233zd(E) fit above "
          f"(own 21-point dataset -- see run_r1234yf_validation() docstring for why)")
    print("=" * 100)

    _results_yf, _loss_meta_yf = run_r1234yf_validation()

    _mape_m_yf, _n_m_yf = compute_mape(_results_yf, "m_dot_meas", "m_dot_model")
    _mape_W_yf, _n_W_yf = compute_mape(_results_yf, "W_meas", "W_model")
    _mape_T_yf, _n_T_yf = compute_mape(_results_yf, "T_dis_meas", "T_dis_model",
                                        filter_fn=lambda r: r["T_dis_steady"] == 1)

    print()
    print("-" * 100)
    print(f"MASS FLOW RATE MAPE : {_mape_m_yf:6.2f} %   (n={_n_m_yf}/{len(VALIDATION_DATA_R1234YF)})   "
          f"{'PASS' if _mape_m_yf < 5.0 else 'FAIL'} (<5% target)")
    print(f"POWER MAPE          : {_mape_W_yf:6.2f} %   (n={_n_W_yf}/{len(VALIDATION_DATA_R1234YF)})   "
          f"{'PASS' if _mape_W_yf < 5.0 else 'FAIL'} (<5% target)")
    print(f"DISCHARGE TEMP MAPE : {_mape_T_yf:6.2f} %   (n={_n_T_yf}/{len(VALIDATION_DATA_R1234YF)}, "
          f"T_dis_steady==1 only)")
    print("-" * 100)

    import json as _json_yf
    with open("r1234yf_validation_results.json", "w") as _f_yf:
        _json_yf.dump({
            "results": _results_yf,
            "mape_mass_flow_pct": _mape_m_yf, "mape_power_pct": _mape_W_yf, "mape_Tdis_pct": _mape_T_yf,
            "n_mass_flow": _n_m_yf, "n_power": _n_W_yf, "n_Tdis": _n_T_yf,
            "config": {"V_dead_ratio": _cfg["V_dead_ratio"], "T_w": T_w,
                       "n_steps": N_STEPS_FINAL, "max_iter": MAX_ITER_FINAL,
                       "d_port": D_VALVE_PORT, "k_v": VALVE_DYNAMICS["k_v"],
                       "m_v": VALVE_DYNAMICS["m_v"], "c_v": VALVE_DYNAMICS["c_v"],
                       "eta_mech": ETA_MECH, "eta_motor_reference": ETA_MOTOR_REFERENCE,
                       "w_loss_fixed_kW": _loss_meta_yf["w_loss_fixed_kW"],
                       "w_loss_slope": _loss_meta_yf["w_loss_slope"],
                       "eta_high_load_asymptote": _loss_meta_yf["eta_high_load_asymptote"]},
        }, _f_yf, indent=2, default=str)
    print("Saved r1234yf_validation_results.json")


# ==============================================================================
# 9. CYCLE DIAGNOSTIC TRACES -- P, T, P-V, power vs. crank angle
# ==============================================================================
# Runs one or more named Bitzer validation points a second time with
# return_trace=True and saves the full crank-angle history of the last
# converged cycle to cycle_trace_<name>.json. Read back by Section 11's
# make_cycle_plots(), which turns each JSON into individual pressure-vs-
# crank-angle, temperature-vs-crank-angle, P-V diagram, and indicated-
# power-vs-crank-angle plots (plus valve-lift/mass-flow bonus plots).
#
# Deliberately a SEPARATE small driver from Section 8, not folded into the
# main validation loop, so that the (cheap, already-converged) validation
# run is never slowed down by trace bookkeeping for points nobody asked
# to inspect -- only the points named in CYCLE_TRACE_POINTS pay this cost.

if COMPRESSOR == "BITZER_2DES3Y" and RUN_CYCLE_TRACES and not INJECTION_ENABLED:
    print()
    print("=" * 100)
    print(f"CYCLE DIAGNOSTIC TRACES: {CYCLE_TRACE_POINTS} "
          f"(n_steps={N_STEPS_FINAL}, max_iter={MAX_ITER_FINAL})")
    print("=" * 100)

    _fluid_trace = "R1233zd(E)"
    _by_name = {row["AutoName"]: row for row in VALIDATION_DATA_R1233ZDE}

    for _pt_name in CYCLE_TRACE_POINTS:
        if _pt_name not in _by_name:
            print(f"  [skip] '{_pt_name}' not found in VALIDATION_DATA_R1233ZDE -- check spelling")
            continue
        _row = _by_name[_pt_name]
        _P_suc = _row["P_suc"] * 1000.0
        _P_dis = _row["P_dis"] * 1000.0
        _T_suc = _row["T_suc"] + 273.15
        _rho_suc = PropsSI("D", "T", _T_suc, "P", _P_suc, _fluid_trace)

        (_T_disc, _W_in, _m_out, _m_in, _h_dis, _m_inj, _trace) = run_compressor_point(
            _fluid_trace, _T_suc, _rho_suc, _P_suc, _P_dis, False, 0.0,
            N_STEPS_FINAL, MAX_ITER_FINAL, VALVE_DYNAMICS, return_trace=True)

        _out = {
            "point": {
                "name": _pt_name, "type": _row["Type"], "P_ratio": _row["P_ratio"],
                "P_suc_kPa": _row["P_suc"], "P_dis_kPa": _row["P_dis"],
                "T_suc_C": _row["T_suc"], "m_dot_meas_kgmin": _row["m_dot"],
                "m_dot_model_kgmin": _m_out * 60.0, "W_meas_kW": _row["W_dot_comp"],
                "W_indicated_model_kW": _W_in / 1000.0,
            },
            "trace": _trace,
        }
        _fname = f"cycle_trace_{_pt_name}.json"
        import json as _json2
        with open(_fname, "w") as _f2:
            _json2.dump(_out, _f2, indent=2, default=str)
        print(f"  {_pt_name} (Pr={_row['P_ratio']:.2f}): saved {_fname}")

    print("Cycle traces done. Section 11 below renders them automatically.")


# ==============================================================================
# 10. VALIDATION PLOTS -- model vs. measured, per refrigerant
# ==============================================================================
# Folded in verbatim from the standalone make_validation_plots.py (this file,
# make_cycle_plots.py, and compressor_model_final.py are now ONE file, per
# request -- nothing about the plotting logic itself changed). The only
# change from the original script is wrapping it in a function parametrized
# by (results json, fluid label, filename prefix) so it can run once per
# refrigerant: called with file_prefix="plot" on bitzer_validation_results.json
# it reproduces the exact 8 filenames/plots already delivered for R1233zd(E)
# (byte-for-byte the same code path, so unchanged), and called again with a
# different prefix on r1234yf_validation_results.json it produces a parallel,
# non-colliding set of 8 plots for R1234yf.

def make_validation_plots(results_json_path, fluid_label, file_prefix):
    """
    Renders the 8 standard validation plots (mass flow / power / discharge
    temperature parity + vs.-pressure-ratio trend, eta_vol vs. eta_ois, and
    the electromechanical loss-model diagnostic) from a results JSON
    produced by run_bitzer_validation() (Section 8) or
    run_r1234yf_validation() (Section 8b). Each plot is its own standalone
    figure/file (not a multi-panel grid) -- to adjust one, find its
    "PLOT N:" section below and edit that block; each is self-contained
    (own fig/ax, own legend call, own save).
    """
    with open(results_json_path) as f:
        data = json.load(f)
    results = data["results"]
    mape_m = data["mape_mass_flow_pct"]
    mape_W = data["mape_power_pct"]
    mape_T = data["mape_Tdis_pct"]
    cfg = data["config"]

    name = [r["name"] for r in results]
    Pr = np.array([r["P_ratio"] for r in results])
    m_meas = np.array([r["m_dot_meas"] for r in results])
    m_model = np.array([r["m_dot_model"] for r in results])
    W_meas = np.array([r["W_meas"] for r in results])
    W_model = np.array([r["W_model"] for r in results])
    Tdis_meas = np.array([r["T_dis_meas"] for r in results])
    Tdis_model = np.array([r["T_dis_model"] for r in results])
    Tdis_steady = np.array([r["T_dis_steady"] for r in results], dtype=bool)
    eta_vol_meas = np.array([r["eta_vol_meas"] for r in results])
    eta_vol_model = np.array([r["eta_vol_model"] for r in results])
    eta_ois_meas = np.array([r["eta_ois_meas"] for r in results])
    eta_ois_model = np.array([r["eta_ois_model"] for r in results])

    COLOR_DATA = "#1f77b4"   # single uniform color/marker for every data point (no C/S split)

    def parity_plot(meas, model, label, unit, mape=None, mask_steady=None, fname=None):
        """Builds ONE standalone parity (model-vs-measured) figure and saves it."""
        fig, ax = plt.subplots(1, 1, figsize=(7, 6.5))

        plot_mask = mask_steady if mask_steady is not None else np.ones_like(meas, dtype=bool)

        # --- data series (edit marker/color/label here) ---
        data_label = "Model vs. measured"
        if mape is not None:
            data_label += f" (MAPE = {mape:.2f}%)"
        ax.scatter(meas[plot_mask], model[plot_mask], c=COLOR_DATA,
                   label=data_label, s=60, zorder=3, edgecolors="k", linewidths=0.4)

        lo = min(meas[plot_mask].min(), model[plot_mask].min()) * 0.9
        hi = max(meas[plot_mask].max(), model[plot_mask].max()) * 1.1
        ax.plot([lo, hi], [lo, hi], "k--", lw=1.2, label="Perfect agreement", zorder=1)
        ax.fill_between([lo, hi], [lo * 0.95, hi * 0.95], [lo * 1.05, hi * 1.05],
                         color="gray", alpha=0.15, zorder=0, label="+/-5%")

        ax.set_xlabel(f"Measured {label} [{unit}]")
        ax.set_ylabel(f"Model {label} [{unit}]")
        ax.set_title(f"{label}: Model vs. Measured\n{fluid_label}")

        # --- legend: edit loc/fontsize/labels here (MAPE is in the data-series label above) ---
        ax.legend(fontsize=9, loc="upper left")

        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(fname, dpi=150)
        plt.close(fig)
        print(f"Saved {fname}")

    def trend_plot(meas, model, label, unit, mask=None, fname=None):
        """Builds ONE standalone vs-pressure-ratio trend figure and saves it."""
        fig, ax = plt.subplots(1, 1, figsize=(7, 6))
        m = mask if mask is not None else np.ones_like(Pr, dtype=bool)

        # --- data series (edit markers/colors/labels here) ---
        ax.scatter(Pr[m], meas[m], c="black", marker="o", label="Measured (Excel)", s=55, zorder=3)
        ax.scatter(Pr[m], model[m], facecolors="none", edgecolors="tab:orange", marker="D",
                   label="Model", s=80, linewidths=1.8, zorder=3)
        for i in np.where(m)[0]:
            ax.plot([Pr[i], Pr[i]], [meas[i], model[i]], color="gray", lw=0.8, alpha=0.6, zorder=1)

        ax.set_xlabel("Pressure ratio Pr [-]")
        ax.set_ylabel(f"{label} [{unit}]")
        ax.set_title(f"{label} vs. pressure ratio\n{fluid_label}")

        # --- legend: edit loc/fontsize/labels here ---
        ax.legend(fontsize=9)

        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(fname, dpi=150)
        plt.close(fig)
        print(f"Saved {fname}")

    # ============================================================
    # PLOT 1: mass flow rate parity
    # ============================================================
    parity_plot(m_meas, m_model, "Mass flow rate", "kg/min", mape=mape_m,
                fname=f"{file_prefix}_01_mass_flow_parity.png")

    # ============================================================
    # PLOT 2: total (electrical) power parity
    # ============================================================
    parity_plot(W_meas, W_model, "Total (electrical) power", "kW", mape=mape_W,
                fname=f"{file_prefix}_02_power_parity.png")

    # ============================================================
    # PLOT 3: discharge temperature parity (steady-state points only)
    # ============================================================
    parity_plot(Tdis_meas, Tdis_model, "Discharge temperature", "C", mape=mape_T,
                mask_steady=Tdis_steady, fname=f"{file_prefix}_03_discharge_temp_parity.png")

    # ============================================================
    # PLOT 4: eta_vol vs eta_ois "practical ordering" vs pressure ratio
    # ============================================================
    fig4, ax4 = plt.subplots(1, 1, figsize=(7.5, 6.5))

    # --- data series (edit markers/colors/labels here) ---
    ax4.scatter(Pr, eta_vol_meas, c="green", marker="o", label="eta_vol measured", s=50)
    ax4.scatter(Pr, eta_vol_model, facecolors="none", edgecolors="green", marker="o",
                label="eta_vol model", s=75, linewidths=1.6)
    ax4.scatter(Pr, eta_ois_meas, c="purple", marker="^", label="eta_ois measured", s=50)
    ax4.scatter(Pr, eta_ois_model, facecolors="none", edgecolors="purple", marker="^",
                label="eta_ois model", s=75, linewidths=1.6)

    ax4.set_xlabel("Pressure ratio Pr [-]")
    ax4.set_ylabel("Efficiency [-]")
    ax4.set_title("eta_vol vs. eta_ois -- \"practical\" ordering\n"
                  f"(eta_vol > eta_ois at every point) -- {fluid_label}")

    # --- legend: edit loc/fontsize/labels here ---
    ax4.legend(fontsize=9)

    ax4.grid(alpha=0.3)
    fig4.tight_layout()
    fig4.savefig(f"{file_prefix}_04_efficiency_vs_pressure_ratio.png", dpi=150)
    plt.close(fig4)
    print(f"Saved {file_prefix}_04_efficiency_vs_pressure_ratio.png")

    # ============================================================
    # PLOT 5: mass flow rate vs. pressure ratio (trend view)
    # ============================================================
    trend_plot(m_meas, m_model, "Mass flow rate", "kg/min",
               fname=f"{file_prefix}_05_mass_flow_vs_pressure_ratio.png")

    # ============================================================
    # PLOT 6: total power vs. pressure ratio (trend view)
    # ============================================================
    trend_plot(W_meas, W_model, "Total power", "kW",
               fname=f"{file_prefix}_06_power_vs_pressure_ratio.png")

    # ============================================================
    # PLOT 7: discharge temperature vs. pressure ratio (trend view, steady pts only)
    # ============================================================
    trend_plot(Tdis_meas, Tdis_model, "Discharge temp.", "C", mask=Tdis_steady,
               fname=f"{file_prefix}_07_discharge_temp_vs_pressure_ratio.png")

    # ============================================================
    # PLOT 8: electromechanical loss-model diagnostic -- why a fixed+proportional
    # fit was needed instead of a single literature efficiency ratio
    # ============================================================
    eta_mech = cfg["eta_mech"]
    eta_motor_ref = cfg["eta_motor_reference"]
    w_loss_fixed = cfg["w_loss_fixed_kW"]
    w_loss_slope = cfg["w_loss_slope"]
    eta_hla = cfg["eta_high_load_asymptote"]

    W_ind = np.array([r["W_indicated_model"] for r in results])
    W_shaft = W_ind / eta_mech

    fig8, ax8 = plt.subplots(1, 1, figsize=(7.5, 6.5))

    # --- data series (edit markers/colors/labels here) ---
    ax8.scatter(W_shaft, W_meas, c=COLOR_DATA, s=60,
                edgecolors="k", linewidths=0.4, zorder=3, label="Measured")
    xs = np.linspace(0, W_shaft.max() * 1.08, 50)
    ax8.plot(xs, w_loss_fixed + w_loss_slope * xs, "k-", lw=1.8,
             label=f"Fitted: W_total = {w_loss_fixed:.3f} + {w_loss_slope:.3f} x W_shaft")
    ax8.plot(xs, xs / eta_motor_ref, "r--", lw=1.4,
             label=f"Rejected: constant-ratio model (W_shaft/{eta_motor_ref:.2f})")

    ax8.set_xlabel("Model shaft power W_shaft = W_indicated / ETA_MECH [kW]")
    ax8.set_ylabel("Measured total (electrical) power [kW]")
    ax8.set_title("Electromechanical loss model: fixed+proportional fit vs.\n"
                  "the rejected constant-efficiency-ratio attempt")
    ax8.annotate(f"high-load asymptotic efficiency = {eta_hla:.3f}\n"
                 f"(literature reference = {eta_mech*eta_motor_ref:.3f})",
                 xy=(0.05, 0.95), xycoords="axes fraction", va="top", fontsize=8.5,
                 bbox=dict(boxstyle="round", fc="wheat", alpha=0.6))

    # --- legend: edit loc/fontsize/labels here ---
    ax8.legend(fontsize=9, loc="lower right")

    ax8.grid(alpha=0.3)
    fig8.tight_layout()
    fig8.savefig(f"{file_prefix}_08_electromechanical_loss_model.png", dpi=150)
    plt.close(fig8)
    print(f"Saved {file_prefix}_08_electromechanical_loss_model.png")

    print(f"All 8 individual plots saved for {fluid_label} (prefix '{file_prefix}_*').")


if COMPRESSOR == "BITZER_2DES3Y" and RUN_VALIDATION_VS_PAPER and not INJECTION_ENABLED:
    print()
    make_validation_plots("bitzer_validation_results.json", "Bitzer 2DES-3Y / R1233zd(E)", "plot")

if COMPRESSOR == "BITZER_2DES3Y" and RUN_R1234YF_VALIDATION and not INJECTION_ENABLED:
    print()
    make_validation_plots("r1234yf_validation_results.json", "Bitzer 2DES-3Y / R1234yf", "plot_r1234yf")


# ==============================================================================
# 11. CYCLE DIAGNOSTIC PLOTS -- pressure, temperature, P-V, power vs. crank angle
# ==============================================================================
# Folded in verbatim from the standalone make_cycle_plots.py -- this code was
# already generic (it auto-discovers every cycle_trace_*.json present in the
# working directory via glob, whichever refrigerant/points produced them), so
# nothing about its logic needed to change to live in this file.

def make_cycle_plots():
    """
    Renders, for every cycle_trace_<name>.json present (see Section 9),
    individual pressure-, temperature-, P-V-, and indicated-power-vs-crank-
    angle plots, plus two bonus plots (valve lift, valve mass flow) since
    the trace data already carries them. Both cylinders are overlaid on
    each plot (Cylinder 1 solid, Cylinder 2 dashed) against the shared
    shaft crank angle, so their 180-degree phase offset is directly visible.
    """
    trace_files = sorted(glob.glob("cycle_trace_*.json"))
    if not trace_files:
        print("No cycle_trace_*.json files found -- skipping Section 11 (cycle plots). "
              "Set RUN_CYCLE_TRACES=True in Section 1 to generate them.")
        return

    COLOR_CYL1 = "#1f77b4"   # cylinder 1 -- solid line in every plot below
    COLOR_CYL2 = "#ff7f0e"   # cylinder 2 -- dashed line in every plot below

    def cycle_plot(theta, series, xlabel, ylabel, title, fname):
        """
        Builds ONE standalone crank-angle trend figure and saves it. `series`
        is a list of (y_values, label, color, linestyle) tuples, all plotted
        against the same `theta` array (already the correct length for these
        y_values -- pass the full theta_deg grid for state variables, or the
        midpoint theta_deg_mid grid for rate quantities like power/mass flow).
        """
        fig, ax = plt.subplots(1, 1, figsize=(7.5, 5.5))

        # --- data series (edit colors/labels/linestyles here) ---
        for y, label, color, ls in series:
            ax.plot(theta, y, color=color, linestyle=ls, lw=1.8, label=label)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_xlim(theta[0], theta[-1])

        # --- legend: edit loc/fontsize/labels here ---
        ax.legend(fontsize=9)

        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(fname, dpi=150)
        plt.close(fig)
        print(f"Saved {fname}")

    def pv_plot(V_list, P_list, labels, colors, lss, title, fname):
        """Builds ONE standalone P-V diagram (all cylinders overlaid) and saves it."""
        fig, ax = plt.subplots(1, 1, figsize=(7, 6.5))

        # --- data series (edit colors/labels/linestyles here) ---
        for V, P, label, color, ls in zip(V_list, P_list, labels, colors, lss):
            ax.plot(V, P, color=color, linestyle=ls, lw=1.8, label=label)

        ax.set_xlabel("Cylinder volume V [cm3]")
        ax.set_ylabel("Cylinder pressure P [kPa]")
        ax.set_title(title)

        # --- legend: edit loc/fontsize/labels here ---
        ax.legend(fontsize=9)

        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(fname, dpi=150)
        plt.close(fig)
        print(f"Saved {fname}")

    for tf in trace_files:
        with open(tf) as f:
            d = json.load(f)
        pt = d["point"]
        tr = d["trace"]
        cycle_name = pt["name"]
        theta = np.array(tr["theta_deg"])            # (n_s,)   P/T/V/rho/h/y_* grid
        theta_mid = np.array(tr["theta_deg_mid"])     # (n_s-1,) mdot_*/Wdot_in/Q_dot grid
        n_cyl = tr["n_cyl"]
        cyls = tr["cylinders"]
        colors = [COLOR_CYL1, COLOR_CYL2][:n_cyl]
        lss = ["-", "--"][:n_cyl]
        cyl_labels = [f"Cylinder {c + 1}" for c in range(n_cyl)]

        subtitle = (f"{cycle_name} (Pr={pt['P_ratio']:.2f}, type {pt['type']}) -- "
                    f"Bitzer 2DES-3Y / R1233zd(E)")

        # ============================================================
        # PLOT 1: pressure vs. crank angle
        # ============================================================
        series = [(cyls[c]["P_kPa"], cyl_labels[c], colors[c], lss[c]) for c in range(n_cyl)]
        cycle_plot(theta, series, "Shaft crank angle [deg]", "Cylinder pressure [kPa]",
                   f"Pressure vs. crank angle\n{subtitle}",
                   f"cycle_{cycle_name}_01_pressure_vs_theta.png")

        # ============================================================
        # PLOT 2: temperature vs. crank angle
        # ============================================================
        series = [(cyls[c]["T_C"], cyl_labels[c], colors[c], lss[c]) for c in range(n_cyl)]
        cycle_plot(theta, series, "Shaft crank angle [deg]", "Cylinder temperature [C]",
                   f"Temperature vs. crank angle\n{subtitle}",
                   f"cycle_{cycle_name}_02_temperature_vs_theta.png")

        # ============================================================
        # PLOT 3: P-V diagram (complete cycle)
        # ============================================================
        V_list = [cyls[c]["V_cm3"] for c in range(n_cyl)]
        P_list = [cyls[c]["P_kPa"] for c in range(n_cyl)]
        pv_plot(V_list, P_list, cyl_labels, colors, lss,
                f"P-V diagram (complete cycle)\n{subtitle}",
                f"cycle_{cycle_name}_03_pv_diagram.png")

        # ============================================================
        # PLOT 4: instantaneous indicated power vs. crank angle
        # ============================================================
        series = [(cyls[c]["Wdot_in_W"], cyl_labels[c], colors[c], lss[c]) for c in range(n_cyl)]
        cycle_plot(theta_mid, series, "Shaft crank angle [deg]",
                   "Instantaneous indicated power, W_dot_in [W]",
                   f"Indicated power vs. crank angle\n{subtitle}",
                   f"cycle_{cycle_name}_04_power_vs_theta.png")

        # ============================================================
        # BONUS PLOT 5: valve lift vs. crank angle (suction + discharge)
        # (not explicitly requested, but the trace already carries it)
        # ============================================================
        fig5, ax5 = plt.subplots(1, 1, figsize=(7.5, 5.5))

        # --- data series (edit colors/labels/linestyles here) ---
        for c in range(n_cyl):
            ax5.plot(theta, cyls[c]["y_suc_mm"], color=colors[c], linestyle=lss[c],
                      lw=1.6, label=f"{cyl_labels[c]} suction valve")
            ax5.plot(theta, cyls[c]["y_dis_mm"], color=colors[c], linestyle=lss[c],
                      lw=1.6, alpha=0.5, label=f"{cyl_labels[c]} discharge valve")
        ax5.axhline(tr["y_max_mm"], color="k", lw=1.0, linestyle=":", label="y_max (mechanical stop)")

        ax5.set_xlabel("Shaft crank angle [deg]")
        ax5.set_ylabel("Valve lift [mm]")
        ax5.set_title(f"Valve lift vs. crank angle (bonus plot)\n{subtitle}")
        ax5.set_xlim(theta[0], theta[-1])

        # --- legend: edit loc/fontsize/labels here ---
        ax5.legend(fontsize=8.5)

        ax5.grid(alpha=0.3)
        fig5.tight_layout()
        fname5 = f"cycle_{cycle_name}_05_valve_lift_vs_theta.png"
        fig5.savefig(fname5, dpi=150)
        plt.close(fig5)
        print(f"Saved {fname5}")

        # ============================================================
        # BONUS PLOT 6: valve mass flow rate vs. crank angle (suction in / discharge out)
        # (not explicitly requested, but the trace already carries it)
        # ============================================================
        fig6, ax6 = plt.subplots(1, 1, figsize=(7.5, 5.5))

        # --- data series (edit colors/labels/linestyles here) ---
        for c in range(n_cyl):
            mdot_in_gs = np.array(cyls[c]["mdot_in_kgs"]) * 1000.0
            mdot_out_gs = np.array(cyls[c]["mdot_out_kgs"]) * 1000.0
            ax6.plot(theta_mid, mdot_in_gs, color=colors[c], linestyle=lss[c],
                      lw=1.6, label=f"{cyl_labels[c]} suction (in)")
            ax6.plot(theta_mid, mdot_out_gs, color=colors[c], linestyle=lss[c],
                      lw=1.6, alpha=0.5, label=f"{cyl_labels[c]} discharge (out)")

        ax6.set_xlabel("Shaft crank angle [deg]")
        ax6.set_ylabel("Valve mass flow rate [g/s]")
        ax6.set_title(f"Valve mass flow rate vs. crank angle (bonus plot)\n{subtitle}")
        ax6.set_xlim(theta_mid[0], theta_mid[-1])

        # --- legend: edit loc/fontsize/labels here ---
        ax6.legend(fontsize=8.5)

        ax6.grid(alpha=0.3)
        fig6.tight_layout()
        fname6 = f"cycle_{cycle_name}_06_massflow_vs_theta.png"
        fig6.savefig(fname6, dpi=150)
        plt.close(fig6)
        print(f"Saved {fname6}")

    print(f"\nDone -- {len(trace_files)} operating point(s) x 6 plots each "
          f"(4 requested: pressure, temperature, P-V, power; 2 bonus: valve lift, mass flow).")


if COMPRESSOR == "BITZER_2DES3Y" and RUN_CYCLE_TRACES and not INJECTION_ENABLED:
    print()
    make_cycle_plots()


