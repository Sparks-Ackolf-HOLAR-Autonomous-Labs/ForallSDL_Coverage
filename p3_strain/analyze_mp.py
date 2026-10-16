#!/usr/bin/env python3
"""Turn one in.mp_kappa run (log + profile) into a kappa JSON point.

Usage: analyze_mp.py <run.log> <run.profile> <out.json>

Same math as pilot/src/mp_kappa.py: kappa = j / <|dT/dz|>, with
j = E_swapped / (2 A t_prod) (two flux directions), slopes fit on both
halves of the running-averaged profile excluding the swap slabs (bins 0
and NBIN/2) and their neighbors. Pure stdlib on purpose: CHPC's default
python3 has no numpy.
"""
import json
import sys

NBIN = 20
EV_PS_A_TO_W_MK = 1.602176634e-19 / (1e-12 * 1e-10)  # eV/(ps A K) -> W/(m K)


def slope(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return sxy / sxx


def main():
    logf, proff, outf = sys.argv[1], sys.argv[2], sys.argv[3]
    meta = {}
    for line in open(logf):
        p = line.split()
        if not p:
            continue
        if p[0] == "MPT_EWARM":
            meta["e_warm"] = float(p[1])
        elif p[0] == "MPT_EPROD":
            meta["e_prod"] = float(p[1])
        elif p[0] == "MPT_BOX":
            meta["box"] = [float(v) for v in p[1:4]]
        elif p[0] == "MPT_NATOMS":
            meta["natoms"] = int(p[1])
        elif p[0] == "MPT_META":
            kv = p[1:]
            meta.update({kv[i]: kv[i + 1] for i in range(0, len(kv), 2)})

    rows = []
    for line in open(proff):
        p = line.split()
        if len(p) == 4 and not line.startswith("#"):
            rows.append((float(p[1]), float(p[3])))
    prof = rows[-NBIN:]  # last running-averaged block: (reduced z, T)
    lx, ly, lz = meta["box"]
    z = [r[0] * lz for r in prof]
    tt = [r[1] for r in prof]
    h = NBIN // 2
    s1 = slope(z[2:h - 1], tt[2:h - 1])
    s2 = slope(z[h + 2:NBIN - 1], tt[h + 2:NBIN - 1])
    grad = 0.5 * (abs(s1) + abs(s2))  # K per Angstrom
    prod_ps = int(meta["prod_steps"]) * 0.001
    de = meta["e_prod"] - meta["e_warm"]  # eV over production
    j = de / (2 * lx * ly * prod_ps)      # eV / (ps A^2)
    kappa = j / grad * EV_PS_A_TO_W_MK

    result = {
        "x": float(meta["x"]), "T": float(meta["T"]),
        "ncz": int(meta["ncz"]), "seed": int(meta["seed"]),
        "method": "Muller-Plathe rNEMD",
        "party": "CHPC-MD", "machine": "notchpeak-guest",
        "natoms": meta["natoms"],
        "box_A": [round(v, 2) for v in meta["box"]],
        "kappa_W_mK": round(kappa, 3),
        "slopes_K_per_A": [round(s1, 5), round(s2, 5)],
        "deltaT_across_half_box_K": round(grad * lz / 2, 1),
        "swapped_eV": round(de, 1),
        "equil_warm_prod_steps": [50000, 100000, 500000],
        # protocol variants (in.mp_kappa2); absent in the original
        # campaign's logs, where these were fixed at 4 and 50
        "ncxy": int(meta.get("ncxy", 4)),
        "swap_every": int(meta.get("swap", 50)),
    }
    with open(outf, "w") as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result))


if __name__ == "__main__":
    main()
