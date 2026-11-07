import json
import os
import subprocess
import sys


from neo2_util import golden_replay


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




def test_libneo_is_the_pin_or_main_of_the_commit_date(tmp_path):
    libneo = new_repo(tmp_path / "libneo")
    early = commit(libneo, {"a": "1"}, "early", "2025-01-01T00:00:00+00:00")
    commit(libneo, {"a": "2"}, "late", "2025-06-01T00:00:00+00:00")
    commit(libneo, {"a": "3"}, "later", "2026-01-01T00:00:00+00:00")
    neo2 = new_repo(tmp_path / "neo2")
    unpinned = commit(neo2, {"x": "1"}, "unpinned", "2025-03-01T00:00:00+00:00")
    commit(neo2, {"cmake/dependency_versions.json": json.dumps({"libneo": early})},
                    "pin", "2025-07-01T00:00:00+00:00")
    assert golden_replay.libneo_revision(neo2, unpinned, libneo) == early
    assert golden_replay.libneo_revision(neo2, "HEAD", libneo) == early


FAKE_SOLVER = """#!{python}
import os, h5py
memory = bytearray({size} * 2**20)
with h5py.File("out.h5", "w") as f:
    f["result"] = [1.0, float(os.environ.get("FAKE_RESULT", "2"))]
"""

CASE_MAKEFILE = """all: summary.h5
summary.h5: .run
\ttest -f shared.txt
\tcp out.h5 summary.h5
.run: neo_2.x
\t./neo_2.x
\ttouch .run
neo_2.x: $(CURRENT_NEO2_QL)
\tln -fs $< $@
clean:
\trm -f .run out.h5 summary.h5 neo_2.x
"""


def toy_case(data, roles=("current",), makefile=CASE_MAKEFILE):
    # Like the real cases, a run directory links a file above the case.
    data.mkdir(parents=True, exist_ok=True)
    (data / "shared.txt").write_text("shared")
    for role in roles:
        (data / "toy" / role).mkdir(parents=True)
        (data / "toy" / role / "Makefile").write_text(makefile)
        (data / "toy" / role / "shared.txt").symlink_to("../../shared.txt")


def fake_tree(path, size):
    for variant in ("ql", "par"):
        binary = golden_replay.executable(path, variant)
        binary.parent.mkdir(parents=True)
        binary.write_text(FAKE_SOLVER.format(python=sys.executable, size=size))
        binary.chmod(0o755)
    (path / "build" / "golden-history.json").write_text(json.dumps(
        {"neo2": "a" * 40, "libneo": "b" * 40, "date": "2026-11-01T00:00:00+00:00",
         "subject": "fake", "build": "ok"}))
    return path


def test_record_measures_paired_runs(tmp_path):
    data = tmp_path / "golden"
    toy_case(data)
    current = fake_tree(tmp_path / "current", 300)
    reference = fake_tree(tmp_path / "reference", 1)
    out = tmp_path / "record"
    result = golden_replay.record(data, current, out, ["toy"], paired=reference, concurrent=True)
    toy = result["cases"]["toy"]
    assert toy["status"] == "ok"
    assert toy["wall_s"] > 0 and toy["reference_wall_s"] > 0
    assert toy["peak_rss_mib"][0] > 300 > toy["reference_peak_rss_mib"][0]
    assert (out / "toy" / "summary.h5").exists()
    assert golden_replay.compare(out, out / "reference", "main") == {"toy": {"/result": 0.0}}
    stored = json.loads((out / "record.json").read_text())
    assert stored["changes"]["main"]["toy"]["/result"] == 0.0


def test_failed_case_is_marked(tmp_path):
    data = tmp_path / "golden"
    toy_case(data, makefile=CASE_MAKEFILE.replace("./neo_2.x", "false"))
    result = golden_replay.record(data, fake_tree(tmp_path / "tree", 1), tmp_path / "out", ["toy"])
    assert result["cases"]["toy"]["status"] == "failed"
