"""Wall time and peak memory of the NEO-2 runs in a golden-record job.

  python -m neo2_util.golden_measure wrap LOG BINARY
      print an executable that runs BINARY through `rank` and logs to LOG
  python -m neo2_util.golden_measure rank LOG BINARY ARGS...
      run one MPI rank or serial stage
  python -m neo2_util.golden_measure summary LOG...
      one row per case and log: runs, median rank-0 wall time, peak memory

The case Makefiles run CURRENT_NEO2_* and REFERENCE_NEO2_*, so pointing these
at the wrapper measures every run without changing a case.
"""
import json
import os
import resource
import signal
import statistics
import subprocess
import sys
import time
from pathlib import Path


def rank(log, binary, args):
    start = time.monotonic()
    child = subprocess.Popen([binary, *args])
    for sig in (signal.SIGTERM, signal.SIGINT):
        signal.signal(sig, lambda s, _: child.send_signal(s))
    code = child.wait()
    peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    # ru_maxrss is in bytes on macOS and in KiB on Linux.
    peak_mib = peak / 2**20 if sys.platform == "darwin" else peak / 2**10
    number = next((os.environ[k] for k in ("OMPI_COMM_WORLD_RANK", "PMI_RANK", "PMIX_RANK")
                   if k in os.environ), "0")
    with open(log, "a") as out:
        out.write(json.dumps({"case": Path.cwd().parent.name, "rank": int(number),
                              "wall_s": round(time.monotonic() - start, 3),
                              "peak_rss_mib": round(peak_mib, 1), "exit": code}) + "\n")
    return code


def wrap(log, binary):
    log, binary = Path(log).resolve(), Path(binary).resolve()
    path = log.parent / f"{log.stem}-{binary.name}"
    path.write_text(f'#!/bin/sh\nexec "{sys.executable}" -m neo2_util.golden_measure '
                    f'rank "{log}" "{binary}" "$@"\n')
    path.chmod(0o755)
    return path


def measurements(log):
    """Per case: summed rank-0 wall time, peak memory per rank and exit status."""
    cases = {}
    lines = Path(log).read_text().splitlines() if Path(log).exists() else []
    for entry in map(json.loads, lines):
        case = cases.setdefault(entry["case"], {"wall_s": 0.0, "peaks": {}, "exit": 0})
        if entry["rank"] == 0:
            case["wall_s"] = round(case["wall_s"] + entry["wall_s"], 3)
        case["peaks"][entry["rank"]] = max(case["peaks"].get(entry["rank"], 0.0),
                                           entry["peak_rss_mib"])
        case["exit"] = case["exit"] or entry["exit"]
    return {name: {"wall_s": c["wall_s"], "exit": c["exit"],
                   "peak_rss_mib": [c["peaks"][r] for r in sorted(c["peaks"])]}
            for name, c in cases.items()}


def summary(logs):
    print("case\tbuild\truns\tmedian_wall_s\tpeak_rss_mib")
    for log in logs:
        entries = [json.loads(line) for line in Path(log).read_text().splitlines()]
        for case in sorted({e["case"] for e in entries}):
            mine = [e for e in entries if e["case"] == case]
            walls = [e["wall_s"] for e in mine if e["rank"] == 0]
            print(f"{case}\t{Path(log).stem}\t{len(walls)}\t"
                  f"{statistics.median(walls) if walls else float('nan'):.1f}\t"
                  f"{max(e['peak_rss_mib'] for e in mine):.0f}")


def main(argv):
    if argv[0] == "wrap":
        print(wrap(argv[1], argv[2]))
    elif argv[0] == "rank":
        return rank(argv[1], argv[2], argv[3:])
    elif argv[0] == "summary":
        summary(argv[1:])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
