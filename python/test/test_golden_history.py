import json
import os
import subprocess

from neo2_util import golden_history, golden_measure, golden_replay
from test_golden_replay import fake_tree, toy_case


def test_collect_measures_runs_of_the_golden_job(tmp_path):
    data = tmp_path / "golden"
    toy_case(data, ("current", "reference"))
    current = fake_tree(tmp_path / "current", 1)
    reference = fake_tree(tmp_path / "reference", 1)
    runs = tmp_path / "runs"
    runs.mkdir()
    for role, log, tree, result in (("current", "current", current, "2"),
                                    ("reference", "base", reference, "4")):
        wrapper = golden_measure.wrap(runs / f"{log}.ranks", golden_replay.executable(tree, "ql"))
        env = dict(os.environ, CURRENT_NEO2_QL=str(wrapper), FAKE_RESULT=result)
        subprocess.run(["make", "-s", "-C", str(data / "toy" / role)], env=env, check=True)
    out = tmp_path / "record.json"
    golden_history.collect(data, current, reference, runs, out)
    record = json.loads(out.read_text())
    toy = record["cases"]["toy"]
    assert toy["status"] == "ok" and toy["reference_wall_s"] > 0
    changes = record["changes"]
    assert changes["previous"]["toy"]["/result"] == changes["main"]["toy"]["/result"] == 0.5


def test_render_lists_events_and_draws_the_panels(tmp_path):
    record = {"neo2": "a" * 40, "libneo": "b" * 40, "date": "2026-11-01T00:00:00+00:00",
              "subject": "fake", "build": "ok",
              "cases": {"toy": {"status": "ok", "wall_s": 2.0, "reference_wall_s": 1.0,
                                "peak_rss_mib": [300.0, 200.0]}},
              "changes": {"main": {"toy": {"/result": 1e-9}}},
              "event": {"source": "neo2", "pr": 1, "description": "toy change"}}
    path = tmp_path / "record.json"
    path.write_text(json.dumps(record))
    golden_history.render([path], tmp_path / "site")
    page = (tmp_path / "site" / "index.html").read_text()
    assert "toy change" in page and "pull/1" in page
    for name in ("changes.png", "changes_main.png", "wall_time.png", "memory.png"):
        assert (tmp_path / "site" / name).stat().st_size > 0
