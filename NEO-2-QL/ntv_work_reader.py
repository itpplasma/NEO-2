"""Independent, strict reader of NEO2_COMPLEX_WORK_V1 streams."""
from pathlib import Path
import math
import struct


def binary_header(path, basis=False):
    path = Path(path)
    size = 24 if basis else 28
    with path.open('rb') as stream:
        data = stream.read(size)
    if len(data) != size:
        raise ValueError('Truncated work header')
    expected = b'N2BAS001' if basis else b'N2WORK01'
    if data[:8] != expected:
        raise ValueError('Wrong work magic/version')
    if data[8:12] == b'\x04\x03\x02\x01':
        endian = '<'
    elif data[8:12] == b'\x01\x02\x03\x04':
        endian = '>'
    else:
        raise ValueError('Invalid native byte-order marker')
    if basis:
        record_size, count = struct.unpack(endian + 'Iq', data[12:])
        if record_size != 164:
            raise ValueError('Wrong pitch-reconstruction record width')
    else:
        real_bits, complex_bits, count = struct.unpack(endian + 'IIq', data[12:])
        if (real_bits, complex_bits) != (64, 128):
            raise ValueError('Unsupported work floating-point representation')
        record_size = 16
    if count < 0 or path.stat().st_size != size + record_size * count:
        raise ValueError('Incomplete or excess work payload')
    return endian, count, size


def complex_values(path):
    endian, count, offset = binary_header(path)
    with Path(path).open('rb') as stream:
        stream.seek(offset)
        for _ in range(count):
            real, imag = struct.unpack(endian + 'dd', stream.read(16))
            if not (math.isfinite(real) and math.isfinite(imag)):
                raise ValueError('Nonfinite complex work payload')
            yield complex(real, imag)


def basis_values(path):
    endian, count, offset = binary_header(path, basis=True)
    with Path(path).open('rb') as stream:
        stream.seek(offset)
        for _ in range(count):
            values = struct.unpack(endian + 'iii19d', stream.read(164))
            if not all(math.isfinite(x) for x in values[3:]):
                raise ValueError('Nonfinite pitch reconstruction')
            yield values


def read_work(path):
    path = Path(path)
    rows = [line.split() for line in path.read_text().splitlines()]
    if not rows or rows[0] != ['NEO2_COMPLEX_WORK_V1']:
        raise ValueError('Wrong work schema')
    if rows[-1] != ['END_NEO2_COMPLEX_WORK_V1']:
        raise ValueError('Incomplete work session')
    metadata = {}
    for row in rows[1:-1]:
        if not row:
            raise ValueError('Empty work metadata row')
        metadata.setdefault(row[0], []).append(row[1:])
    length = int(metadata['VECTOR_LENGTH'][0][0])
    vectors = {}
    for role, count, filename in metadata['VECTOR']:
        if role in vectors:
            raise ValueError('Duplicate work-vector role')
        vector_path = path.parent / filename
        _, actual, _ = binary_header(vector_path)
        if actual != length or int(count) != length:
            raise ValueError('Work-vector layout mismatch')
        vectors[role] = vector_path
    basis = path.parent / metadata['BASIS_FILE'][0][0]
    binary_header(basis, basis=True)
    return dict(metadata=metadata, vectors=vectors, basis=basis, length=length)
