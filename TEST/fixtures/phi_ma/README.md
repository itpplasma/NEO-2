Synthetic circular QL smoke fixture, with no machine or experiment data.

Provenance: the hand-written `smoke.py` from the neo2-phi1 W1 inspection
(`/mnt/storage/neo2-phi1/w1-20261022/smoke.py`, 2026-10-22). The background has
R=3 m, a=0.3 m, iota=0.7 and B=2 T*(1-r/R*cos(theta)), at four synthetic
surfaces s=0.1,0.3,0.6,0.9 with r=a*sqrt(s). The perturbation file declares
n=1 with zero magnetic coefficients. The species, profiles and rotation are
artificial. A2 is zero to exercise the scalar J_hat units boundary.

This coarse fixture checks the production source and current plumbing and
linear symmetries. It is not a converged cylindrical susceptibility oracle.
