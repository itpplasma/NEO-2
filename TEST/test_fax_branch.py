"""Physical branch labels survive either sign of the toroidal field."""
import importlib.util
from pathlib import Path
import sys
import numpy as np

spec = importlib.util.spec_from_file_location('fax', Path(sys.argv[1]))
fax = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fax)
d = dict(lag=0,num_species=1,npassing=np.array([1]),offset=np.array([0]),
         eta=np.array([0.,0.4]),bhat=np.array([1.]),
         g_real=np.arange(4.,dtype=float).reshape(4,1,1)*np.ones((1,3,1)),
         phi_direction=np.array(1.))
x, y = fax.pitch_rows(d,0)
d['phi_direction'] = np.array(-1.)
u, v = fax.pitch_rows(d,0)
np.testing.assert_array_equal(u,-x[::-1])
np.testing.assert_array_equal(v,y[:,::-1])
d.update(bhat=np.array([1.,1.5,2.,1.5,1.]),
         dlogb_dphi=np.array([0.,1.,0.,-1.,0.]))
np.testing.assert_allclose(fax.position(d),[0.,-np.pi/2,-np.pi,np.pi/2,0.],atol=1e-15)
print('PASS: both physical parallel-branch orientations')
