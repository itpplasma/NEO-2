"""Relative output changes between golden-record runs, and the history of
golden-record refreshes and tolerance changes.

    python -m neo2_util.golden_changes events DATA_REPO NEO2_REPO
"""

import json
import re
import subprocess
import sys
from pathlib import Path

GOLDEN_PATH = "TESTS/NEO-2/golden_record"
# Run metadata that differs between any two runs, as in the cases' compare_h5.py.
METADATA = {
    "CMake_Blas_Lib", "CMake_Build_Type", "CMake_Flags_Debug",
    "CMake_Flags_Release", "CMake_SuiteSparse_Dir", "CMake_System",
    "MPILib Version", "NEO-2 Version", "NEO-2 Version Additional", "Run_0",
    "Timestamp_start_0",
}
# Multistage runs add one run number and start and stop time per stage.
STAGE_METADATA = re.compile(r"^(Run|Timestamp_start|Timestamp_stop)_[0-9]+$")
TOLERANCE = re.compile(r"(?i)\b(\w*tol|\w*accuracy|delta_\w+)\s*\??[:=]\s*([0-9][0-9.e+-]*)")


def git(repo, *args):
    return subprocess.run(["git", "-C", str(repo), *args], check=True,
                          capture_output=True, text=True).stdout.strip()


def relative_changes(reference, other):
    """max|other - reference| / max|reference| for every numeric dataset.

    This is the measure of the golden-record comparison (neo2_util.hdf5tools).
    A dataset with an absolute 'accuracy' attribute counts as unchanged
    within that accuracy, as in the comparison. A dataset present in only one
    file, or of different shape, has no relative change and is stored as None.
    """
    import h5py
    import numpy
    from .hdf5tools import _comparison_keys

    changes = {}

    def walk(a, b):
        ka, kb = _comparison_keys(a), _comparison_keys(b)
        for key in sorted(set(ka) | set(kb)):
            if key in METADATA or STAGE_METADATA.match(key):
                continue
            name = f"{a.name.rstrip('/')}/{key}"
            if key not in ka or key not in kb:
                changes[name] = None
                continue
            x, y = a[ka[key]], b[kb[key]]
            if isinstance(x, h5py.Group) and isinstance(y, h5py.Group):
                walk(x, y)
            elif isinstance(x, h5py.Group) or isinstance(y, h5py.Group) or x.shape != y.shape:
                changes[name] = None
            elif x.dtype.kind in "biufc" and y.dtype.kind in "biufc":
                u, v = numpy.asarray(x[()], dtype=complex), numpy.asarray(y[()], dtype=complex)
                if u.size == 0:
                    continue
                scale = numpy.max(numpy.abs(u)) or 1.0
                delta = numpy.max(numpy.abs(v - u))
                accuracy = x.attrs.get("accuracy", -1.0)
                if numpy.isfinite(accuracy) and 0 <= accuracy and delta <= accuracy:
                    delta = 0.0
                diff = float(delta / scale)
                changes[name] = diff if numpy.isfinite(diff) else 1.0
            else:
                changes[name] = 0.0 if numpy.array_equal(x[()], y[()]) else 1.0

    with h5py.File(reference, "r") as a, h5py.File(other, "r") as b:
        walk(a, b)
    return changes



def events(data_repo, neo2_repo):
    """Commits that refresh golden records or change a comparison tolerance."""
    rows = []
    for repo, source, paths in ((data_repo, "data", [GOLDEN_PATH]),
                                (neo2_repo, "neo2", ["ci", ".github/workflows"])):
        log = git(repo, "log", "--reverse", "--format=%x00%H%x09%cs%x09%s", "--name-only", "--", *paths)
        for block in filter(None, log.split("\0")):
            head, *files = block.strip("\n").split("\n")
            sha, date, subject = head.split("\t", 2)
            files = [f for f in files if f]
            records = sorted({Path(f).relative_to(GOLDEN_PATH).parts[0] for f in files
                              if source == "data" and re.search(r"/(legacy|reference)/summary\.h5$", f)})
            patch = git(repo, "show", "--format=", "-U0", sha, "--",
                        *[f for f in files if not f.endswith(".h5")]) if files else ""
            tolerance = _tolerances(patch, "-") != _tolerances(patch, "+")
            kind = "+".join(k for k, hit in (("record", records), ("tolerance", tolerance)) if hit)
            if kind:
                rows.append({"date": date, "source": source, "commit": sha, "kind": kind,
                             "cases": ",".join(records), "subject": subject})
    return sorted(rows, key=lambda r: r["date"])


def _tolerances(patch, sign):
    # Moved or reformatted settings cancel; only changed values count.
    lines = [l[1:] for l in patch.splitlines() if l.startswith(sign) and not l.startswith(sign * 3)]
    return sorted(m for l in lines for m in TOLERANCE.findall(l))




def main(argv=None):
    args = sys.argv[1:] if argv is None else argv
    if len(args) == 3 and args[0] == "events":
        print(json.dumps(events(args[1], args[2]), indent=1))
        return 0
    print(__doc__, file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
