#!/usr/bin/env python3
"""Results, run time and memory of NEO-2 at every golden-record change.

One command builds a NEO-2/libneo pair and runs the golden-record cases:

    python -m neo2_util.golden_history run --neo2 REV --data GOLDEN --out OUT

OUT/record.json holds run times, peak memory per MPI rank and relative
changes; OUT/<case>/summary.h5 holds the absolute outputs, which stay with the
golden-record data because some cases use non-public inputs. The other
subcommands build, record, compare, list events and render the site.
"""

import argparse
import datetime
import json
import os
import re
import resource
import shutil
import signal
import statistics
import subprocess
import sys
import threading
import time
from pathlib import Path

CASES = ("lorentz", "ql", "par", "performance_ql", "transport_k")
LIBNEO_URL = "https://github.com/itpplasma/libneo.git"
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


def libneo_revision(neo2_repo, neo2_rev, libneo_repo):
    """The pinned libneo, or before the pin the libneo main of the commit date."""
    try:
        pins = json.loads(git(neo2_repo, "show", f"{neo2_rev}:cmake/dependency_versions.json"))
        return git(libneo_repo, "rev-parse", pins["libneo"] + "^{commit}")
    except subprocess.CalledProcessError:
        pass
    date = git(neo2_repo, "log", "-1", "--format=%cI", neo2_rev)
    rev = git(libneo_repo, "rev-list", "-1", "--first-parent", f"--before={date}",
              "refs/remotes/origin/main" if _has_ref(libneo_repo, "refs/remotes/origin/main")
              else "main")
    if not rev:
        raise ValueError(f"libneo has no main commit before {date}")
    return rev


def _has_ref(repo, ref):
    return subprocess.run(["git", "-C", str(repo), "show-ref", "--verify", "--quiet", ref]).returncode == 0


def build(source, neo2_rev, libneo_source, libneo_rev, work, jobs=None):
    """Check out the pair under work/ and build Release; return its source root."""
    work = Path(work).resolve()
    neo2 = work / "neo2"
    libneo = work / "libneo"
    for path, repo, rev in ((neo2, source, neo2_rev), (libneo, libneo_source, libneo_rev)):
        if not path.exists():
            subprocess.run(["git", "clone", "-q", "--no-checkout", str(repo), str(path)], check=True)
        git(path, "checkout", "-q", "--detach", rev)
    info = _revisions(neo2, libneo)
    info["build"] = "failed"
    # FETCHCONTENT_SOURCE_DIR_LIBNEO overrides the libneo fetch of every
    # NEO-2 revision since its FetchContent build, pinned or not.
    configure = ["cmake", "-S", str(neo2), "-B", str(neo2 / "build"), "-G", "Ninja",
                 "-DCMAKE_BUILD_TYPE=Release", f"-DFETCHCONTENT_SOURCE_DIR_LIBNEO={libneo}"]
    compile_ = ["cmake", "--build", str(neo2 / "build")] + (["-j", str(jobs)] if jobs else [])
    with open(work / "build.log", "w") as log:
        ok = all(subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT).returncode == 0
                 for cmd in (configure, compile_))
    if ok and all(executable(neo2, v).exists() for v in ("ql", "par")):
        info["build"] = "ok"
    (neo2 / "build" / "golden-history.json").parent.mkdir(parents=True, exist_ok=True)
    (neo2 / "build" / "golden-history.json").write_text(json.dumps(info, indent=1))
    return neo2


def executable(source, variant):
    return Path(source) / "build" / f"NEO-2-{variant.upper()}" / f"neo_2_{variant}.x"


def describe(source):
    """Revisions of a built tree, from its build note or from its checkouts."""
    note = Path(source) / "build" / "golden-history.json"
    if note.exists():
        return json.loads(note.read_text())
    info = _revisions(source, Path(source) / "build" / "_deps" / "libneo-src")
    info["build"] = "ok"
    return info


def _revisions(neo2, libneo):
    # A pair exists from the later of its two commit dates on.
    info = {"neo2": git(neo2, "rev-parse", "HEAD"), "neo2_date": git(neo2, "log", "-1", "--format=%cI"),
            "subject": git(neo2, "log", "-1", "--format=%s"), "libneo": "", "libneo_date": ""}
    if Path(libneo).exists():
        info["libneo"] = git(libneo, "rev-parse", "HEAD")
        info["libneo_date"] = git(libneo, "log", "-1", "--format=%cI")
    info["date"] = max(filter(None, (info["neo2_date"], info["libneo_date"])),
                       key=datetime.datetime.fromisoformat)
    return info


def rank(log, binary, args):
    """Run one MPI rank or serial stage and append its time and peak memory."""
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
        out.write(json.dumps({"cwd": os.getcwd(), "rank": int(number),
                              "peak_rss_mib": round(peak_mib, 1),
                              "wall_s": round(time.monotonic() - start, 3), "exit": code}) + "\n")
    return code


def wrap(log, binary):
    """Write an executable that runs binary through rank(); return its path.

    The case Makefiles link neo_2.x to CURRENT_NEO2_* or REFERENCE_NEO2_*, so
    pointing these at the wrapper measures every run without changing a case.
    """
    log, binary = Path(log).resolve(), Path(binary).resolve()
    path = log.parent / f"{log.stem}-{binary.name}"
    path.write_text(f'#!/bin/sh\nexec "{sys.executable}" -m neo2_util.golden_history rank '
                    f'"{log}" "{binary}" "$@"\n')
    path.chmod(0o755)
    return path


def measurements(log):
    """Per case (the directory above the run directory): solver time and peaks.

    The time is the sum over solver runs of rank 0, so it leaves out the
    Python collection steps of a case.
    """
    cases = {}
    lines = Path(log).read_text().splitlines() if Path(log).exists() else []
    for entry in map(json.loads, lines):
        case = cases.setdefault(Path(entry["cwd"]).parent.name,
                                {"wall_s": 0.0, "peaks": {}, "exit": 0})
        if entry["rank"] == 0:
            case["wall_s"] = round(case["wall_s"] + entry["wall_s"], 3)
        case["peaks"][entry["rank"]] = max(case["peaks"].get(entry["rank"], 0.0),
                                           entry["peak_rss_mib"])
        case["exit"] = case["exit"] or entry["exit"]
    return {name: {"wall_s": c["wall_s"], "exit": c["exit"],
                   "peak_rss_mib": [c["peaks"][r] for r in sorted(c["peaks"])]}
            for name, c in cases.items()}


def run_case(data, case, source, work):
    """Run one case of a built tree in a fresh copy; return its measurements."""
    work = Path(work).resolve()
    case_dir = work / case
    shutil.rmtree(case_dir, ignore_errors=True)
    # Cases link files above their directory (../../collect_h5.py), so copy
    # the link targets rather than the links.
    shutil.copytree(Path(data) / case, case_dir)
    log = work / f"{case}.ranks"
    log.unlink(missing_ok=True)
    env = dict(os.environ)
    for variant in ("ql", "par"):
        env[f"CURRENT_NEO2_{variant.upper()}"] = str(wrap(log, executable(source, variant)))
    current = case_dir / "current"
    with open(work / f"{case}.log", "w") as out:
        subprocess.run(["make", "-C", str(current), "clean"], env=env, stdout=out,
                       stderr=subprocess.STDOUT)
        code = subprocess.run(["make", "-C", str(current)], env=env, stdout=out,
                              stderr=subprocess.STDOUT).returncode
    entry = measurements(log).get(case, {"wall_s": 0.0, "exit": 1, "peak_rss_mib": []})
    ok = code == 0 and entry.pop("exit") == 0 and (current / "summary.h5").exists()
    entry.pop("exit", None)
    return {"status": "ok" if ok else "failed", **entry}, current / "summary.h5"


def collect(data, source, reference, history, out):
    """Public record of a golden-record run whose solvers ran through wrap().

    The current runs log to history/current.ranks, the reference runs of the
    same job to history/reference.ranks, and the case directories in data
    hold the outputs of both.
    """
    result = describe(source)
    result["reference"] = describe(reference)["neo2"]
    result["cases"] = {}
    result["changes"] = {"previous": {}, "main": {}}
    current = measurements(Path(history) / "current.ranks")
    paired = measurements(Path(history) / "reference.ranks")
    for case, entry in sorted(current.items()):
        summary = Path(data) / case / "current" / "summary.h5"
        entry["status"] = "ok" if entry.pop("exit") == 0 and summary.exists() else "failed"
        if case in paired:
            entry["reference_wall_s"] = paired[case]["wall_s"]
            entry["reference_peak_rss_mib"] = paired[case]["peak_rss_mib"]
            base = Path(data) / case / "reference" / "summary.h5"
            if summary.exists() and base.exists():
                # The reference is main just before the change, so it is the
                # previous event and main at once.
                changes = relative_changes(base, summary)
                result["changes"]["previous"][case] = result["changes"]["main"][case] = changes
        result["cases"][case] = entry
    Path(out).write_text(json.dumps(result, indent=1))
    return result


def record(data, source, out, cases=CASES, paired=None, concurrent=False):
    """Run the cases of one build, optionally paired with a reference build."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    result = describe(source)
    result["cases"] = {}
    result["changes"] = {}
    if result["build"] != "ok":
        (out / "record.json").write_text(json.dumps(result, indent=1))
        return result
    if paired:
        result["reference"] = describe(paired)["neo2"]
    for case in cases:
        runs = {}

        def go(label, tree):
            runs[label] = run_case(data, case, tree, out / "work" / label)

        jobs = [("current", source)] + ([("reference", paired)] if paired else [])
        if concurrent:
            threads = [threading.Thread(target=go, args=job) for job in jobs]
            for thread in threads:
                thread.start()
            for thread in threads:
                thread.join()
        else:
            for job in jobs:
                go(*job)
        entry, summary = runs["current"]
        if "reference" in runs:
            entry["reference_wall_s"] = runs["reference"][0]["wall_s"]
            entry["reference_peak_rss_mib"] = runs["reference"][0]["peak_rss_mib"]
        if summary.exists():
            (out / case).mkdir(exist_ok=True)
            shutil.copy2(summary, out / case / "summary.h5")
        if "reference" in runs and runs["reference"][1].exists():
            (out / "reference" / case).mkdir(parents=True, exist_ok=True)
            shutil.copy2(runs["reference"][1], out / "reference" / case / "summary.h5")
        result["cases"][case] = entry
    # Keep the solver logs of failed runs; the copied cases are rebuilt at will.
    for log in (out / "work").glob("*/*.log"):
        (out / "logs").mkdir(exist_ok=True)
        shutil.copy2(log, out / "logs" / f"{log.parent.name}-{log.name}")
    shutil.rmtree(out / "work", ignore_errors=True)
    (out / "record.json").write_text(json.dumps(result, indent=1))
    return result


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


def compare(record_dir, other_dir, label):
    """Store the relative changes of record_dir against other_dir under label."""
    path = Path(record_dir) / "record.json"
    result = json.loads(path.read_text())
    changes = {}
    for case in result.get("cases", {}):
        mine, theirs = Path(record_dir) / case / "summary.h5", Path(other_dir) / case / "summary.h5"
        if mine.exists() and theirs.exists():
            changes[case] = relative_changes(theirs, mine)
    result.setdefault("changes", {})[label] = changes
    path.write_text(json.dumps(result, indent=1))
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


def render(records, out):
    """Write the public site: relative changes, run time and memory over time."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rows = sorted((json.loads(Path(p).read_text()) for p in records), key=lambda r: r["date"])
    dates = {id(r): datetime.datetime.fromisoformat(r["date"]) for r in rows}
    cases = sorted({c for r in rows for c in r.get("cases", {})})
    # Paired times: ratio to the reference run on the same node, times the
    # median reference time of the case, which removes drifts in node speed.
    ref = {c: statistics.median([r["cases"][c]["reference_wall_s"] for r in rows
                                 if c in r.get("cases", {}) and "reference_wall_s" in r["cases"][c]]
                                or [1.0]) for c in cases}

    def change(case, entry, r, label):
        values = r.get("changes", {}).get(label, {}).get(case)
        numbers = [v for v in (values or {}).values() if v is not None]
        return max(max(numbers, default=0.0), 1e-17) if values else None

    def wall(case, entry, r):
        if entry.get("reference_wall_s"):
            return entry["wall_s"] / entry["reference_wall_s"] * ref[case]
        return entry["wall_s"]

    panels = [
        ("changes.png", "Largest relative change against the previous event", "relative change", True,
         lambda c, e, r: change(c, e, r, "previous")),
        ("changes_main.png", "Largest relative change against current main", "relative change", True,
         lambda c, e, r: change(c, e, r, "main")),
        ("wall_time.png", "Wall time, paired against a reference on the same node", "wall time (s)",
         False, wall),
        ("memory.png", "Peak resident memory of the largest MPI rank", "peak memory (MiB)", False,
         lambda c, e, r: max(e["peak_rss_mib"], default=None)),
    ]
    for name, title, ylabel, log, value in panels:
        fig, ax = plt.subplots(figsize=(9, 4))
        for case in cases:
            points = [(dates[id(r)], value(case, r["cases"][case], r), r) for r in rows
                      if r.get("cases", {}).get(case, {}).get("status") == "ok"]
            points = [p for p in points if p[1] is not None]
            if not points:
                continue
            line, = ax.plot([p[0] for p in points], [p[1] for p in points], "-", lw=0.8, label=case)
            for when, y, r in points:
                source = r.get("event", {}).get("source")
                if source == "libneo":
                    ax.plot(when, y, "^", color=line.get_color(), mfc="none")
                elif source:
                    ax.plot(when, y, "o", color=line.get_color())
        for r in rows:
            if r["build"] != "ok" or any(c.get("status") != "ok" for c in r["cases"].values()):
                ax.axvline(dates[id(r)], color="0.6", ls=":", lw=0.8)
        if log:
            ax.set_yscale("log")
        ax.set_title(title)
        ax.set_xlabel("commit date")
        ax.set_ylabel(ylabel)
        if ax.has_data():
            ax.legend(fontsize="small", ncol=len(cases))
        fig.autofmt_xdate()
        fig.tight_layout()
        fig.savefig(out / name, dpi=120)
        plt.close(fig)
    _write_html(rows, out)


def _write_html(rows, out):
    import html
    repo = "https://github.com/itpplasma/NEO-2"
    body = []
    for r in reversed(rows):
        event = r.get("event", {})
        failed = [c for c, e in r["cases"].items() if e.get("status") != "ok"]
        status = "build failed" if r["build"] != "ok" else (
            "failed: " + ", ".join(failed) if failed else "ok")
        pr = event.get("pr")
        issue = event.get("issue")
        cells = [r["date"][:10], html.escape(event.get("source", "")),
                 f'<a href="{repo}/commit/{r["neo2"]}">{r["neo2"][:8]}</a>'
                 + (f' (<a href="{repo}/pull/{pr}">#{pr}</a>)' if pr else ""),
                 f'<a href="https://github.com/itpplasma/libneo/commit/{r["libneo"]}">{r["libneo"][:8]}</a>'
                 if r.get("libneo") else "",
                 html.escape(event.get("description", r.get("subject", ""))),
                 html.escape(event.get("classification", ""))
                 + (f' (<a href="{repo}/issues/{issue}">#{issue}</a>)' if issue else ""),
                 status]
        body.append("<tr>" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>")
    figures = "".join(f'<img src="{n}" alt="{n[:-4]}">' for n in
                      ("changes.png", "changes_main.png", "wall_time.png", "memory.png"))
    (out / "index.html").write_text(f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>NEO-2 golden-record history</title>
<style>body{{font-family:sans-serif;max-width:60rem;margin:auto;padding:0 1rem}}
img{{max-width:100%}}table{{border-collapse:collapse;font-size:90%}}
td,th{{border-bottom:1px solid #ccc;padding:.2rem .4rem;text-align:left;vertical-align:top}}</style>
</head><body><h1>NEO-2 golden-record history</h1>
<p>Each row is a build of NEO-2 with the libneo it used, run on the golden-record test cases.
Relative change is max|x &minus; x<sub>ref</sub>| / max|x<sub>ref</sub>| over each output dataset,
largest over the datasets of a case. Filled circles mark NEO-2 events, open triangles libneo events,
dotted lines builds or cases that failed.</p>
{figures}
<table><tr><th>date</th><th>event</th><th>NEO-2</th><th>libneo</th><th>change</th>
<th>classification</th><th>status</th></tr>
{''.join(body)}
</table></body></html>
""")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)

    def case_options(p):
        p.add_argument("--data", required=True, help="golden_record directory of plasma/data")
        p.add_argument("--out", required=True)
        p.add_argument("--cases", nargs="+", default=list(CASES))
        p.add_argument("--paired-build", help="reference NEO-2 tree run on the same node")
        p.add_argument("--concurrent", action="store_true",
                       help="start the build and its reference at the same moment")

    p = sub.add_parser("run", help="build a NEO-2/libneo pair and record it")
    p.add_argument("--neo2", required=True)
    p.add_argument("--libneo", help="default: the pin, or libneo main at the commit date")
    p.add_argument("--source", default=".")
    p.add_argument("--libneo-source", default=LIBNEO_URL)
    p.add_argument("--work", required=True)
    p.add_argument("--jobs", type=int)
    case_options(p)
    p = sub.add_parser("build", help="build a NEO-2/libneo pair")
    for name, kw in (("--neo2", {"required": True}), ("--libneo", {}), ("--source", {"default": "."}),
                     ("--libneo-source", {"default": LIBNEO_URL}), ("--work", {"required": True}),
                     ("--jobs", {"type": int})):
        p.add_argument(name, **kw)
    p = sub.add_parser("record", help="run the cases of a built tree")
    p.add_argument("--build", required=True, help="NEO-2 source tree with build/")
    case_options(p)
    p = sub.add_parser("compare", help="add relative changes against another record")
    p.add_argument("record")
    p.add_argument("other")
    p.add_argument("--label", required=True, choices=("previous", "main"))
    p = sub.add_parser("events", help="list golden-record refreshes and tolerance changes")
    p.add_argument("--data-repo", required=True)
    p.add_argument("--neo2-repo", default=".")
    p = sub.add_parser("render", help="write the public site from record.json files")
    p.add_argument("records", nargs="+")
    p.add_argument("--out", required=True)
    p = sub.add_parser("wrap", help="print a measuring wrapper for a NEO-2 executable")
    p.add_argument("log")
    p.add_argument("binary")
    p = sub.add_parser("collect", help="public record of a golden-record run with wrappers")
    p.add_argument("--data", required=True)
    p.add_argument("--build", required=True)
    p.add_argument("--reference-build", required=True)
    p.add_argument("--history", required=True, help="directory of current.ranks, reference.ranks")
    p.add_argument("--out", required=True)
    p = sub.add_parser("rank")
    p.add_argument("log")
    p.add_argument("binary")
    p.add_argument("args", nargs=argparse.REMAINDER)
    a = parser.parse_args(argv)

    if a.command == "rank":
        return rank(a.log, a.binary, a.args)
    if a.command == "wrap":
        print(wrap(a.log, a.binary))
        return 0
    if a.command == "collect":
        collect(a.data, a.build, a.reference_build, a.history, a.out)
        return 0
    if a.command in ("run", "build"):
        libneo_source = a.libneo_source
        work = Path(a.work)
        if not Path(libneo_source).exists():
            mirror = work / "libneo"
            if not mirror.exists():
                subprocess.run(["git", "clone", "-q", "--no-checkout", libneo_source, str(mirror)],
                               check=True)
            libneo_source = mirror
        libneo = a.libneo or libneo_revision(a.source, a.neo2, libneo_source)
        tree = build(a.source, a.neo2, libneo_source, libneo, work, a.jobs)
        print(json.dumps(describe(tree)))
        if a.command == "build":
            return 0 if describe(tree)["build"] == "ok" else 1
        a.build = tree
    if a.command in ("run", "record"):
        result = record(a.data, a.build, a.out, a.cases, a.paired_build, a.concurrent)
        print(json.dumps(result.get("cases", {}), indent=1))
        return 0
    if a.command == "compare":
        changes = compare(a.record, a.other, a.label)
        for case, values in changes.items():
            print(case, max((v for v in values.values() if v is not None), default=0.0), sep="\t")
    if a.command == "events":
        columns = ("date", "source", "commit", "kind", "cases", "subject")
        print(*columns, sep="\t")
        for row in events(a.data_repo, a.neo2_repo):
            print(*(row[c] for c in columns), sep="\t")
    if a.command == "render":
        render(a.records, a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
