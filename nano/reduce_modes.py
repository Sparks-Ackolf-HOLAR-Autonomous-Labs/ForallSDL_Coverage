#!/usr/bin/env python3
"""Reduce fetched almaBTE phononinfo files to compact committed evidence.

Input: a staging dir (default ~/nano_stage/modes/x<X>/) holding
<base>_<T>K.phononinfo, *.Tsweep and sha256sums as written on CHPC.
Output (committed): modes/x<X>.npz with float32 omega, |v|^2 (T-independent)
and C, tau per temperature, plus modes/raw_sha256.json mapping every raw
file to the sha256 computed on CHPC and re-checked here. The raw files
(7 MB each) stay on CHPC home and in the staging dir; the hashes let anyone
holding them prove the npz came from them.
"""
import glob, hashlib, json, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STAGE = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/nano_stage/modes")
OUT = os.path.join(HERE, "modes")
os.makedirs(OUT, exist_ok=True)
TS = [300, 500, 700, 900, 1100]


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


hashes = {}
for xd in sorted(glob.glob(os.path.join(STAGE, "x*"))):
    x = os.path.basename(xd)[1:]
    remote = {}
    if os.path.exists(os.path.join(xd, "sha256sums")):
        for line in open(os.path.join(xd, "sha256sums")):
            h, f = line.split()
            remote[os.path.basename(f)] = h
    arrs, omega, v2 = {}, None, None
    for T in TS:
        fs = glob.glob(os.path.join(xd, f"*_{T}K.phononinfo"))
        if len(fs) != 1:
            sys.exit(f"x={x} T={T}: expected one phononinfo, found {fs}")
        f = fs[0]
        h = sha(f)
        if remote.get(os.path.basename(f), h) != h:
            sys.exit(f"hash mismatch {f}")
        hashes[f"x{x}/{os.path.basename(f)}"] = h
        d = np.loadtxt(f, delimiter=",", skiprows=1)
        w, vv = d[:, 5], (d[:, 8:11] ** 2).sum(1)
        if omega is None:
            omega, v2 = w, vv
        elif not (np.allclose(w, omega) and np.allclose(vv, v2)):
            sys.exit(f"x={x}: dispersion differs between temperatures")
        arrs[f"C_{T}"] = d[:, 6].astype(np.float32)
        arrs[f"tau_{T}"] = d[:, 7].astype(np.float32)
    ts = {}
    for f in glob.glob(os.path.join(xd, "*.Tsweep")):
        hashes[f"x{x}/{os.path.basename(f)}"] = sha(f)
        for line in open(f).read().splitlines()[1:]:
            p = line.split(",")
            ts[str(int(float(p[0])))] = {"rta_xx": float(p[1]), "bte_xx": float(p[2])}
    np.savez_compressed(os.path.join(OUT, f"x{x}.npz"), omega=omega.astype(np.float32),
                        v2=v2.astype(np.float32), **arrs)
    json.dump(ts, open(os.path.join(OUT, f"x{x}_tsweep.json"), "w"), indent=1)
    print("reduced x", x, "modes", len(omega), "tsweep T", sorted(ts))
json.dump(hashes, open(os.path.join(OUT, "raw_sha256.json"), "w"), indent=1, sort_keys=True)
