import json
import os
import subprocess
import sys

import h5py
import pytest

from neo2_util import golden_history


def git(repo, *args, date=None):
    # Fixture repositories must not pick up hooks or signing from user config.
    env = dict(os.environ, GIT_CONFIG_GLOBAL=os.devnull, GIT_CONFIG_NOSYSTEM="1", GIT_AUTHOR_NAME="t", GIT_AUTHOR_EMAIL="t@t",
               GIT_COMMITTER_NAME="t", GIT_COMMITTER_EMAIL="t@t")
    if date:
        env.update(GIT_AUTHOR_DATE=date, GIT_COMMITTER_DATE=date)
    return subprocess.run(["git", "-C", str(repo), *args], check=True, env=env,
                          capture_output=True, text=True).stdout.strip()


def commit(repo, files, message, date="2026-01-01T12:00:00+00:00"):
    for name, text in files.items():
        path = repo / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        git(repo, "add", name)
    git(repo, "commit", "-q", "-m", message, date=date)
    return git(repo, "rev-parse", "HEAD")


def new_repo(path):
    path.mkdir()
    git(path, "init", "-q", "-b", "main")
    return path


def write_h5(path, values):
    with h5py.File(path, "w") as f:
        for key, value in values.items():
            f[key] = value


def test_relative_changes_follow_the_golden_comparison(tmp_path):
    write_h5(tmp_path / "a.h5", {"fulltransp/avnabpsi": [1.0, 2.0], "x": [4.0, -2.0],
                                 "neo2_config/metadata/NEO-2 Version": "a", "gone": 1.0})
    write_h5(tmp_path / "b.h5", {"fulltransp/av_nabla_stor": [1.0, 2.0], "x": [4.0, -2.0 + 4e-9],
                                 "neo2_config/metadata/NEO-2 Version": "b", "new": 1.0})
    changes = golden_history.relative_changes(tmp_path / "a.h5", tmp_path / "b.h5")
    assert changes["/x"] == pytest.approx(1e-9)
    assert changes["/fulltransp/av_nabla_stor"] == 0.0
    assert changes["/gone"] == changes["/new"] == 1.0
    assert not any("Version" in name for name in changes)


def test_events_are_record_refreshes_and_tolerance_changes(tmp_path):
    data = new_repo(tmp_path / "data")
    case = golden_history.GOLDEN_PATH + "/ql"
    commit(data, {f"{case}/Makefile": "TOL ?= 1e-15\n", f"{case}/legacy/summary.h5": "1"}, "add")
    commit(data, {f"{case}/Makefile": "TOL ?= 1e-11\n"}, "relax")
    commit(data, {f"{case}/README.md": "text"}, "docs")
    commit(data, {f"{case}/legacy/summary.h5": "2"}, "refresh")
    neo2 = new_repo(tmp_path / "neo2")
    commit(neo2, {"ci/golden.sh": "compare base-ql ql TOL=1e-15\n"}, "golden")
    commit(neo2, {"ci/golden.sh": "compare base-ql ql TOL=1e-12\n"}, "relax ci")
    commit(neo2, {"ci/run.sh": "echo\n"}, "other")
    rows = golden_history.events(data, neo2)
    found = [(r["source"], r["subject"], r["kind"], r["cases"]) for r in rows]
    assert found == [("data", "add", "record+tolerance", "ql"), ("data", "relax", "tolerance", ""),
                     ("data", "refresh", "record", "ql"), ("neo2", "golden", "tolerance", ""),
                     ("neo2", "relax ci", "tolerance", "")]


def test_libneo_is_the_pin_or_main_of_the_commit_date(tmp_path):
    libneo = new_repo(tmp_path / "libneo")
    early = commit(libneo, {"a": "1"}, "early", "2025-01-01T00:00:00+00:00")
    commit(libneo, {"a": "2"}, "late", "2025-06-01T00:00:00+00:00")
    commit(libneo, {"a": "3"}, "later", "2026-01-01T00:00:00+00:00")
    neo2 = new_repo(tmp_path / "neo2")
    unpinned = commit(neo2, {"x": "1"}, "unpinned", "2025-03-01T00:00:00+00:00")
    commit(neo2, {"cmake/dependency_versions.json": json.dumps({"libneo": early})},
                    "pin", "2025-07-01T00:00:00+00:00")
    assert golden_history.libneo_revision(neo2, unpinned, libneo) == early
    assert golden_history.libneo_revision(neo2, "HEAD", libneo) == early


FAKE_SOLVER = """#!{python}
import os, h5py
memory = bytearray({size} * 2**20)
with h5py.File("out.h5", "w") as f:
    f["result"] = [1.0, float(os.environ.get("FAKE_RESULT", "2"))]
"""

CASE_MAKEFILE = """all: summary.h5
summary.h5: .run
\tcp out.h5 summary.h5
.run: neo_2.x
\t./neo_2.x
\ttouch .run
neo_2.x: $(CURRENT_NEO2_QL)
\tln -fs $< $@
clean:
\trm -f .run out.h5 summary.h5 neo_2.x
"""


def fake_tree(path, size):
    for variant in ("ql", "par"):
        binary = golden_history.executable(path, variant)
        binary.parent.mkdir(parents=True)
        binary.write_text(FAKE_SOLVER.format(python=sys.executable, size=size))
        binary.chmod(0o755)
    (path / "build" / "golden-history.json").write_text(json.dumps(
        {"neo2": "a" * 40, "libneo": "b" * 40, "date": "2026-11-01T00:00:00+00:00",
         "subject": "fake", "build": "ok"}))
    return path


def test_record_measures_paired_runs_and_renders(tmp_path):
    data = tmp_path / "golden"
    (data / "toy" / "current").mkdir(parents=True)
    (data / "toy" / "current" / "Makefile").write_text(CASE_MAKEFILE)
    current = fake_tree(tmp_path / "current", 300)
    reference = fake_tree(tmp_path / "reference", 1)
    out = tmp_path / "record"
    result = golden_history.record(data, current, out, ["toy"], paired=reference, concurrent=True)
    toy = result["cases"]["toy"]
    assert toy["status"] == "ok"
    assert toy["wall_s"] > 0 and toy["reference_wall_s"] > 0
    assert toy["peak_rss_mib"][0] > 300 > toy["reference_peak_rss_mib"][0]
    assert (out / "toy" / "summary.h5").exists()
    assert golden_history.compare(out, out / "reference", "main") == {"toy": {"/result": 0.0}}
    stored = json.loads((out / "record.json").read_text())
    assert stored["changes"]["main"]["toy"]["/result"] == 0.0
    stored["event"] = {"source": "neo2", "pr": 1, "description": "toy change"}
    (out / "record.json").write_text(json.dumps(stored))
    golden_history.render([out / "record.json"], tmp_path / "site")
    page = (tmp_path / "site" / "index.html").read_text()
    assert "toy change" in page and "pull/1" in page
    assert (tmp_path / "site" / "memory.png").stat().st_size > 0


def test_failed_case_is_marked(tmp_path):
    data = tmp_path / "golden"
    (data / "toy" / "current").mkdir(parents=True)
    (data / "toy" / "current" / "Makefile").write_text(CASE_MAKEFILE.replace("./neo_2.x", "false"))
    result = golden_history.record(data, fake_tree(tmp_path / "tree", 1), tmp_path / "out", ["toy"])
    assert result["cases"]["toy"]["status"] == "failed"


def test_collect_measures_runs_of_the_golden_job(tmp_path):
    data = tmp_path / "golden"
    for role in ("current", "reference"):
        (data / "toy" / role).mkdir(parents=True)
        (data / "toy" / role / "Makefile").write_text(CASE_MAKEFILE)
    current = fake_tree(tmp_path / "current", 1)
    reference = fake_tree(tmp_path / "reference", 1)
    history = tmp_path / "history"
    history.mkdir()
    for role, tree, result in (("current", current, "2"), ("reference", reference, "4")):
        wrapper = golden_history.wrap(history / f"{role}.ranks", golden_history.executable(tree, "ql"))
        env = dict(os.environ, CURRENT_NEO2_QL=str(wrapper), FAKE_RESULT=result)
        subprocess.run(["make", "-s", "-C", str(data / "toy" / role)], env=env, check=True)
    out = tmp_path / "record.json"
    golden_history.collect(data, current, reference, history, out)
    toy = json.loads(out.read_text())["cases"]["toy"]
    assert toy["status"] == "ok" and toy["reference_wall_s"] > 0
    changes = json.loads(out.read_text())["changes"]
    assert changes["previous"]["toy"]["/result"] == changes["main"]["toy"]["/result"] == 0.5
