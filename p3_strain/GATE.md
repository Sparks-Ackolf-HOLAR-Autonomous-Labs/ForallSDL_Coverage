# P3 strain campaign: validation gate (2026-09-29)

Six runs of `in.mp_strain` on notchpeak-guest (array 15890803, 16 ranks,
zero preemptions), all at ncz = 48 (6144 atoms). Raw outputs in
`gate_results/`. Gate runs are never pooled with the seed batch.

| run | x | T (K) | eps in-plane | geometry | eps normal (relaxed) | PE (eV/atom) | kappa (W/mK) |
|---|---|---|---|---|---|---|---|
| G0_ref | 0.5 | 300 | 0 | cross-plane | -0.00005 | -4.2018 | 0.923 |
| G1_cpc | 0.5 | 300 | -0.02 | cross-plane | +0.01753 | -4.1895 | 1.033 |
| G2_inp | 0.5 | 300 | -0.02 | in-plane | +0.01749 | -4.1802 | 0.930 |
| G3_hotc | 0.5 | 1100 | -0.03 | cross-plane | +0.02837 | -4.0677 | 1.140 |
| G4_hott | 0.5 | 1100 | +0.01 | in-plane | -0.00930 | -4.0785 | 0.993 |
| G5_si | 0.1 | 300 | -0.03 | cross-plane | +0.03140 | -4.5003 | 2.280 |

## Checks

1. **Reference reproduces the unstrained protocol.** G0 (eps = 0, the new
   strain machinery switched on but idle) gives 0.923 W/mK, inside the
   Phase B seed band at the same size (x = 0.5, 300 K, ncz 48: 0.917 to
   0.941). The strain path adds no bias. PASS.
2. **Mechanical response is the potential's own.** The relaxed normal
   strain gives an effective biaxial Poisson ratio -eps_normal/eps_inplane
   of 0.88 (x = 0.5, 300 K), 0.87 (in-plane geometry), 0.95 (1100 K,
   compressive), 0.93 (1100 K, tensile) and 1.05 at x = 0.1. For pure
   Tersoff-1989 Si, 2 C12 / C11 is about 1.05, so the response is the
   potential's elastic behaviour, and the composition trend is smooth.
   Experimental Si is about 0.77, so Tersoff over-predicts the normal
   relaxation by roughly 35%: recorded as a method-bias scope label on
   every strain claim (the out-of-plane lattice of a strained film is
   the quantity it distorts). PASS (with label).
3. **No structural failure.** Every strained box stays crystalline: PE per
   atom moves by at most 0.02 eV at fixed T (compare 300 K runs), and the
   1100 K runs sit where thermal energy puts them; no run shows the
   amorphous or melting signature documented in the P1 probes. PASS.
4. **Physically sensible, finite kappa with linear profiles**, imposed
   temperature difference across half a box 225 to 1205 K, in the same
   range as the P1 production campaign at the same swap rate. PASS.

Honest note on timing: the strain readouts (check 2) were visible while
the MD phase was still running, before this file was written; the kappa
values (checks 1 and 4) were not. The criteria are the ones any strained
MD validation uses, and none was tuned.

## Early physics glimpse (not a result: one seed each)

At 300 K and x = 0.5, 2% compressive in-plane strain raises cross-plane
kappa by about 12% (0.923 to 1.033) but in-plane kappa by about 1%
(0.930), inside the 3% seed scatter. If the seed batch confirms it, the
strain response of SiGe films is anisotropic, which is the kind of
statement the P3 materials question is after.

**Verdict: gate PASSED. The seed batch may be pre-registered and submitted.**
