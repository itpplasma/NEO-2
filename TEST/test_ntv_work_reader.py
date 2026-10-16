"""Closed-form complex-work and exact strided-ABI reconstruction oracles."""
from pathlib import Path
import importlib.util
import sys
import struct

reader_path = Path(__file__).resolve().parents[1] / 'NEO-2-QL/ntv_work_reader.py'
spec = importlib.util.spec_from_file_location('work_reader', reader_path)
reader = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reader)
folder = Path(sys.argv[1])
data = reader.read_work(folder / 'NEO2_WORK_s0_j1.txt')
n = data['length']
assert n == 8201
expected = {
    'FLUX_k1_b0': lambda j: complex(j, 2*j-3),
    'RHS_k1_b0': lambda j: complex(3*j+7, j*j-1),
    'RHS_k3_b0': lambda j: complex(3*j+7, j*j-1),
    'SOL_k1_b0': lambda j: complex(-j+2, 4*j+1),
    'SOL_k3_b0': lambda j: complex(2*j+3, -j+5),
}
assert data['vectors'].keys() == expected.keys()
for role, function in expected.items():
    assert list(reader.complex_values(data['vectors'][role])) == [function(j) for j in range(1, n+1)]
sum_j = n*(n+1)//2
sum_j2 = n*(n+1)*(2*n+1)//6
q1 = (7*sum_j2-8*sum_j-3*n)/2
q3 = (16*sum_j-15*n)/2
native = data['metadata']['QFLUX'][0]
assert [float(x) for x in native[1:]] == [q1, q3]
basis = list(reader.basis_values(data['basis']))
assert len(basis) == 2 and [b[:3] for b in basis] == [(0, 1, 1), (1, 2, 1)]
assert basis[0][3:6] == (0.5, 1.0, 2.0)
assert basis[0][6:] == tuple(i+10*j for j in range(1, 5) for i in range(1, 5))
# Strict readers must reject damaged data, not silently accept a partial export.
broken = folder / 'damaged.bin'
original = data['vectors']['FLUX_k1_b0'].read_bytes()
broken.write_bytes(original[:-1])
try:
    reader.binary_header(broken)
except ValueError:
    pass
else:
    raise AssertionError('Truncated vector accepted')
# Independently poison one active complex value; strict content reader must fail.
endian, _, offset = reader.binary_header(data['vectors']['FLUX_k1_b0'])
poisoned = bytearray(original)
poisoned[offset:offset+8] = struct.pack(endian+'d', float('nan'))
broken.write_bytes(poisoned)
try:
    list(reader.complex_values(broken))
except ValueError:
    pass
else:
    raise AssertionError('Nonfinite active value accepted')
broken.unlink()
print('All tests passed! Exact complex/strided work and pitch-basis streams.')
