#!/usr/bin/env python3
"""Write nano_design.json: the pre-registered design of the nanostructured
SiGe + light-doping campaign (P1 system, run 2026-09-29).

Method in one line: almaBTE (first-principles force constants, VCA alloy)
gives every phonon mode's heat capacity C, lifetime tau and group velocity
v via `phononinfo`; extra scattering channels are added mode by mode with
Matthiessen's rule and the lattice conductivity is re-summed in RTA:

    kappa_L = (1/3) sum_modes C |v|^2 / (1/tau + |v|/d + 1/tau_dop(omega))

grain boundaries: gray diffuse Casimir rate |v|/d (d = grain size);
dopants: Tamura point-defect rate (pi/6) Gamma omega^2 Dtilde(omega), with
Dtilde the per-atom phonon DOS (integral 3) built from the same modes and
Gamma = c [(1 - M_d/Mbar)^2 + eps (1 - r_d/rbar)^2] (Abeles form).

Everything below is fixed before any production phononinfo run exists.
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))

design = {
    "campaign": "nano-sige-doping-2026-09-29",
    "paper": "P1",
    "operator_go": "issue 55, 2026-09-29: 'go for Nanostructured SiGe and light doping'",
    "question": (
        "Which grain size and composition window minimises the lattice thermal "
        "conductivity of a SiGe thermoelectric leg, how much germanium can "
        "nanostructuring save, and is light doping (B, P, C) a useful extra "
        "thermal lever on top of alloying and grain boundaries?"),
    "tool": {"almaBTE": "1.3.2 (CHPC module almabte/1.3.2)",
             "force_constants": "almaBTE database Si.tar.xz, Ge.tar.xz (as party 4a)",
             "grid": "24x24x24 (as party 4a)",
             "solver": "RTA per-mode re-summation of phononinfo output"},
    "grid": {
        "x_Ge": [0.0, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.85, 1.0],
        "new_vca_builds": [0.1, 0.3],
        "T_K": [300, 500, 700, 900, 1100],
        "grain_size_nm": "log-spaced 5 nm to 10 um, 25 points, plus bulk (no boundary term)",
        "dopants": {
            "none": [0],
            "B": {"n_cm3": [1e19, 3e19, 1e20, 3e20]},
            "P": {"n_cm3": [1e19, 3e19, 1e20, 3e20]},
            "C": {"at_frac": [0.0025, 0.005, 0.01, 0.02]},
        },
    },
    "physics_constants": {
        "mass_amu": {"Si": 28.0855, "Ge": 72.630, "B": 10.811, "P": 30.974, "C": 12.011},
        "tetrahedral_covalent_radius_A": {"Si": 1.17, "Ge": 1.22, "B": 0.88, "P": 1.10, "C": 0.77},
        "radius_source": "Pauling tetrahedral covalent radii",
        "lattice_A": {"Si": 5.431, "Ge": 5.658, "rule": "Vegard in x"},
        "host_mass_and_radius": "Vegard average of Si and Ge (the VCA host)",
        "strain_parameter_eps": {"lower_bound": 0.0, "upper_bound": 39.0,
                                 "upper_source": "Abeles 1963 fitted value for Si-Ge"},
        "boundary_model": "gray diffuse Casimir: rate |v|/d, specularity 0",
        "dos_broadening": "Gaussian, sigma = 0.5% of max omega",
    },
    "gates_N0": {
        "G1_postprocessor": ("per-mode RTA bulk kappa reproduces kappa_Tsweep RTA "
                             "within 0.5% at every (x, T) where party-4a Tsweep "
                             "output exists (x in {0,0.15,0.2,0.5,0.7,0.85,1}, "
                             "T in {300,700,1100})"),
        "G2_tamura_implementation": ("Tamura post-processor applied to the pure-Si "
                                     "modes at 300 K with the x=0.2 VCA mass variance "
                                     "(g = sum f_i (1 - M_i/Mbar)^2) gives kappa within "
                                     "a factor 1.5 of almaBTE's own VCA x=0.2 kappa "
                                     "(different dispersions, so a sanity bound only)"),
        "G3_boundary_limits": ("kappa(d) strictly decreasing in d-> smaller; "
                               "kappa(10 um) within 30% of bulk for alloys; "
                               "pure Si kappa(5 nm) < 25% of bulk"),
        "on_failure": "record gate_failure on the ledger; no N1 result is reported as a finding",
    },
    "predictions": {
        "Pr1": "x=0.2, 300 K, d=20 nm, undoped: kappa_L in [2, 5] W/mK",
        "Pr2": "nanostructuring suppression kappa(20 nm)/kappa_bulk at 300 K is smaller (stronger) at x=0.2 than at x=0.5",
        "Pr3": ("light doping is a weak lever in the alloy: for x in [0.15, 0.85] and "
                "n <= 1e20 cm^-3, B or P lowers kappa_L by less than 10% at every d "
                "under the mass-only bound (eps=0)"),
        "Pr4": ("boundary and dopant scattering show diminishing returns: synergy "
                "S = kappa_both*kappa_bulk/(kappa_dop*kappa_gb) >= 1 at x=0.2, 300 K, "
                "d=20 nm, B 1e20"),
        "Pr5": ("open-data check (N3): the effective grain size inferred from each "
                "Starrydata sintered/nanostructured Si-rich (x in [0.15,0.25]) 300 K "
                "curve, treating total kappa as an upper bound on kappa_L, has a "
                "median in [10, 300] nm"),
    },
    "deliverables": {
        "D1": "kappa_L(x, d, T) table plus a compact law with declared validity domain",
        "D2": ("Ge-savings frontier: minimum x reaching kappa_L <= 2.5 and <= 4.0 W/mK "
               "at 300 K and 1100 K, as a function of grain size"),
        "D3": "dopant lever table: percent change of kappa_L for B, P, C vs n, x, d, both eps bounds",
        "D4": "screening decision: (x, d >= 10 nm) minimising kappa_L at 300 K and at 1100 K",
    },
    "declared_limits": [
        "VCA plus Tamura mass disorder; no short-range order, no local force-constant change",
        "gray diffuse boundary with one length; real grain size distributions and interface roughness spectra are not modelled",
        "dopants scatter phonons only by mass and size; electron-phonon scattering at 1e20 cm^-3 and dopant precipitation are out of scope",
        "strain term for B, P, C is bracketed (eps 0 to 39), not computed; a MACE local-relaxation gate is the named follow-up to narrow it",
        "RTA; pure Si and Ge RTA underestimates full BTE by about 4-12% (party 4a), alloys by about 1.5%",
        "no electronic thermal conductivity, so comparisons with measured total kappa are one-sided bounds",
    ],
    "deviation_policy": ("infrastructure reruns allowed with a ledger deviation entry; "
                         "any change to gates, grid, constants or predictions after "
                         "registration is a new pre-registration, not an edit"),
}

if __name__ == "__main__":
    with open(os.path.join(HERE, "nano_design.json"), "w") as f:
        json.dump(design, f, indent=1, sort_keys=True)
        f.write("\n")
    print("wrote nano_design.json")
