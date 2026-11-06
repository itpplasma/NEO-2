import json
import os
import subprocess
import sys

import h5py
import pytest

from neo2_util import golden_changes


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
                                 "neo2_config/metadata/NEO-2 Version": "a", "gone": 1.0,
                                 "neo2_config/metadata/Timestamp_stop_2": "a", "dbp": [1e-12, 1.0]})
    with h5py.File(tmp_path / "a.h5", "r+") as f:
        f["dbp"].attrs["accuracy"] = 1e-10
    write_h5(tmp_path / "b.h5", {"fulltransp/av_nabla_stor": [1.0, 2.0], "x": [4.0, -2.0 + 4e-9],
                                 "neo2_config/metadata/NEO-2 Version": "b", "new": 1.0,
                                 "neo2_config/metadata/Timestamp_stop_2": "b", "dbp": [5e-11, 1.0]})
    changes = golden_changes.relative_changes(tmp_path / "a.h5", tmp_path / "b.h5")
    assert changes["/x"] == pytest.approx(1e-9)
    assert changes["/fulltransp/av_nabla_stor"] == 0.0
    assert changes["/gone"] is None and changes["/new"] is None
    assert changes["/dbp"] == 0.0
    assert not any("Version" in name or "Timestamp" in name for name in changes)


def test_events_are_record_refreshes_and_tolerance_changes(tmp_path):
    data = new_repo(tmp_path / "data")
    case = golden_changes.GOLDEN_PATH + "/ql"
    commit(data, {f"{case}/Makefile": "TOL ?= 1e-15\n", f"{case}/legacy/summary.h5": "1"}, "add")
    commit(data, {f"{case}/Makefile": "TOL ?= 1e-11\n"}, "relax")
    commit(data, {f"{case}/README.md": "text"}, "docs")
    commit(data, {f"{case}/compare_h5.py": "f(a, delta_relative=1e-8)\n"}, "compare")
    commit(data, {f"{case}/compare_h5.py": "f(a,  delta_relative=1e-8)\n"}, "reformat")
    commit(data, {f"{case}/compare_h5.py": "f(a, delta_relative=2e-8)\n"}, "relax par")
    commit(data, {f"{case}/legacy/summary.h5": "2"}, "refresh")
    neo2 = new_repo(tmp_path / "neo2")
    commit(neo2, {"ci/golden.sh": "compare base-ql ql TOL=1e-15\n"}, "golden")
    commit(neo2, {"ci/golden.sh": "compare base-ql ql TOL=1e-12\n"}, "relax ci")
    commit(neo2, {"ci/run.sh": "echo\n"}, "other")
    rows = golden_changes.events(data, neo2)
    found = [(r["source"], r["subject"], r["kind"], r["cases"]) for r in rows]
    assert found == [("data", "add", "record+tolerance", "ql"), ("data", "relax", "tolerance", ""),
                     ("data", "compare", "tolerance", ""),
                     ("data", "relax par", "tolerance", ""), ("data", "refresh", "record", "ql"), ("neo2", "golden", "tolerance", ""),
                     ("neo2", "relax ci", "tolerance", "")]


