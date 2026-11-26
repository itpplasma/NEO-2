"""Build past NEO-2/libneo revisions and run the golden-record cases on them.

    python -m neo2_util.golden_replay run --neo2 REV --data GOLDEN --work W --out OUT

OUT/record.json holds run times, peak memory per MPI rank and the relative
changes against other records; OUT/<case>/summary.h5 holds the absolute
outputs, which stay with the golden-record data because some cases use
non-public inputs. A NEO-2 commit is built with libneo at its pin, or before
pins existed with libneo main at the commit date.
"""

import argparse
import datetime
import json
import os
import shutil
import subprocess
import sys
import threading
from pathlib import Path

from .golden_changes import git, relative_changes
from .golden_measure import measurements, wrap

CASES = ("lorentz", "ql", "par", "performance_ql", "transport_k")
LIBNEO_URL = "https://github.com/itpplasma/libneo.git"
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
    failed = entry.pop("exit")
    ok = code == 0 and failed == 0 and (current / "summary.h5").exists()
    return {"status": "ok" if ok else "failed", **entry}, current / "summary.h5"


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

    for name in ("run", "build"):
        p = sub.add_parser(name, help=f"{name} a NEO-2/libneo pair")
        p.add_argument("--neo2", required=True)
        p.add_argument("--libneo", help="default: the pin, or libneo main at the commit date")
        p.add_argument("--source", default=".")
        p.add_argument("--libneo-source", default=LIBNEO_URL)
        p.add_argument("--work", required=True)
        p.add_argument("--jobs", type=int)
        if name == "run":
            case_options(p)
    p = sub.add_parser("record", help="run the cases of a built tree")
    p.add_argument("--build", required=True, help="NEO-2 source tree with build/")
    case_options(p)
    p = sub.add_parser("compare", help="add relative changes against another record")
    p.add_argument("record")
    p.add_argument("other")
    p.add_argument("--label", required=True, choices=("previous", "main"))
    a = parser.parse_args(argv)

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
    changes = compare(a.record, a.other, a.label)
    print(json.dumps(changes, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
