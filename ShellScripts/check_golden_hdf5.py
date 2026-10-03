"""Check golden dataset structure, finiteness, and complex values."""

import argparse

import h5py
import numpy as np


def datasets(file):
    result = {}

    def visit(name, item):
        if isinstance(item, h5py.Dataset):
            result[name] = item[()]

    file.visititems(visit)
    return result


def check(reference, current, tolerance):
    with h5py.File(reference) as file:
        left = datasets(file)
    with h5py.File(current) as file:
        right = datasets(file)
    if left.keys() != right.keys():
        raise ValueError(f"Dataset paths differ: {left.keys() ^ right.keys()}")
    for name, value in left.items():
        other = right[name]
        if np.shape(value) != np.shape(other):
            raise ValueError(f"Dataset shape differs: {name}")
        for array in (value, other):
            if np.issubdtype(np.asarray(array).dtype, np.number):
                if not np.isfinite(array).all():
                    raise ValueError(f"Non-finite dataset: {name}")
        # The existing comparator handles real values, but skips complex data.
        if np.iscomplexobj(value) or np.iscomplexobj(other):
            scale = np.max(np.abs(value)) if np.size(value) else 0.0
            if scale == 0.0:
                scale = 1.0
            if (np.abs(value - other) > tolerance * scale).any():
                raise ValueError(f"Complex dataset differs: {name}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference")
    parser.add_argument("current")
    parser.add_argument("tolerance", type=float)
    args = parser.parse_args()
    check(args.reference, args.current, args.tolerance)
