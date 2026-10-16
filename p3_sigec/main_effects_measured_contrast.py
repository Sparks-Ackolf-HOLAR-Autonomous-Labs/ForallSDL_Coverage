#!/usr/bin/env python3
"""Main-effect share of var(ln kappa_L) over the full 9000-point P3 factorial, registered stiffness contrast vs the measured radius-dependent contrast. Not pre-registered; appended to results/p3_sigec_measured_contrast.json."""
import sys, numpy as np, json, os
sys.path.insert(0, "."); 
import sensitivity_measured_contrast as S
P=S.P
ym=np.log(P.response("measured")); ys=np.log(np.load("results/response_stiffness.npy"))
FULL=np.array(P.FULL)
def me(y):
    v=y.var(); out={}
    for j,n in enumerate(P.NAMES):
        means=[y[FULL[:,j]==l].mean() for l in range(P.NL[j])]
        cnt=[(FULL[:,j]==l).sum() for l in range(P.NL[j])]
        out[n]=float(round(100*np.average((np.array(means)-y.mean())**2,weights=cnt)/v,2))
    return out
r={"main_effect_pct_of_var_lnkappa_full_factorial":{"stiffness_registered":me(ys),"measured_curve":me(ym)}}
print(r)
f="results/p3_sigec_measured_contrast.json"; d=json.load(open(f)); d.update(r); json.dump(d,open(f,"w"),indent=1)
