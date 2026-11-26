"""Golden-record history site: the record of a CI golden change and the plots.

A golden-change PR leaves a record of its merge-queue run (collect); records
of replayed commits come from golden_replay. render draws relative changes,
paired run time and peak memory over the commit date and lists the events.
"""

import argparse
import datetime
import json
import statistics
import sys
from pathlib import Path

from .golden_changes import relative_changes
from .golden_measure import measurements
from .golden_replay import describe


def collect(data, source, reference, runs, out):
    """Public record of a golden-record run whose solvers ran through wrap().

    The current runs log to runs/current.ranks, the runs of the merge base in
    the same job to runs/base.ranks, and the case directories in data hold
    the outputs of both.
    """
    result = describe(source)
    result["reference"] = describe(reference)["neo2"]
    result["cases"] = {}
    result["changes"] = {"previous": {}, "main": {}}
    current = measurements(Path(runs) / "current.ranks")
    paired = measurements(Path(runs) / "base.ranks")
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
    p = sub.add_parser("collect", help="record of a golden-record CI run with measured solvers")
    p.add_argument("--data", required=True, help="golden_record directory the cases ran in")
    p.add_argument("--build", required=True, help="NEO-2 tree of the change")
    p.add_argument("--reference-build", required=True, help="NEO-2 tree of the merge base")
    p.add_argument("--runs", required=True, help="directory of current.ranks and base.ranks")
    p.add_argument("--out", required=True)
    p = sub.add_parser("render", help="write the site from record.json files")
    p.add_argument("records", nargs="+")
    p.add_argument("--out", required=True)
    a = parser.parse_args(argv)
    if a.command == "collect":
        collect(a.data, a.build, a.reference_build, a.runs, a.out)
    else:
        render(a.records, a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
