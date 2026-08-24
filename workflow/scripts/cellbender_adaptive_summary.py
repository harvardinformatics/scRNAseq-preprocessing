#!/usr/bin/env python
"""Aggregate the per-sample CellBender adaptive-run status files into one run-level summary.

Writes a TSV (one row per sample) and prints a short report: how many samples needed a re-run,
which were resolved by halving the learning rate, and which were not (kept the initial run)."""

import argparse
from collections import Counter
from pathlib import Path

COLUMNS = [
    "sample", "outcome", "initial_learning_rate", "initial_verdict",
    "reran", "rerun_learning_rate", "rerun_verdict", "kept_run",
]


def read_status(path):
    record = {}
    for line in Path(path).read_text().splitlines():
        if "\t" in line:
            key, value = line.split("\t", 1)
            record[key] = value
    return record


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", required=True)
    ap.add_argument("status_files", nargs="*")
    args = ap.parse_args(argv)

    rows = sorted((read_status(p) for p in args.status_files),
                  key=lambda r: r.get("sample", ""))

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as handle:
        handle.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            handle.write("\t".join(row.get(col, "NA") for col in COLUMNS) + "\n")

    counts = Counter(row.get("outcome", "NA") for row in rows)
    reran = [row["sample"] for row in rows if row.get("reran") == "true"]
    unresolved = [row["sample"] for row in rows
                  if row.get("outcome") == "RERUN_DID_NOT_RESOLVE"]

    print(f"[cellbender-adaptive-summary] {len(rows)} sample(s):")
    for outcome in ("NORMAL_FIRST_TRY", "RERUN_RESOLVED", "RERUN_DID_NOT_RESOLVE",
                    "RERUN_SUGGESTED_BUT_ADAPTIVE_OFF"):
        if counts.get(outcome):
            print(f"  {outcome}: {counts[outcome]}")
    if reran:
        print(f"  re-ran with halved learning rate: {', '.join(reran)}")
    if unresolved:
        print(f"  re-run did NOT resolve the learning curve (kept initial): "
              f"{', '.join(unresolved)}")
    print(f"[cellbender-adaptive-summary] wrote {args.output}")


if __name__ == "__main__":
    main()
