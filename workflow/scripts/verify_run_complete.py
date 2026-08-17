#!/usr/bin/env python3
"""Verify a workflow run produced its required outputs, excusing low-quality samples.

Snakemake's own process exit code can be unreliable with the SLURM executor: a terminal job
failure (e.g. an OOM that exhausts retries) does not always propagate to a non-zero exit, so
the runner batch job can show COMPLETED when real work failed. This script gives the runner
job a correct final state independent of that.

It checks whether every required target (the rule-all inputs, passed via --targets-file)
exists, and fails ONLY when a required output is missing for a sample that was NOT flagged as
low quality. Failures confined to flagged low-quality samples - which are legitimately
quarantined and cannot be diagnosed without running the workflow - are excused, so a run whose
only problems are low-quality samples still reports success. A missing output for any other
sample (a real failure) fails the run.

Exit status: 0 if the run is complete for every non-low-quality sample; 1 otherwise. Intended
to be invoked from the Snakefile onsuccess/onerror handlers, whose sys.exit() propagates this
status to the runner job.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


def read_lines(path: Path) -> list[str]:
    if not path.is_file():
        return []
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def read_sample_ids(samplesheet: Path) -> list[str]:
    if not samplesheet.is_file():
        return []
    ids: list[str] = []
    with samplesheet.open() as fh:
        next(fh, None)  # header
        for line in fh:
            if line.strip():
                ids.append(line.split("\t")[0].strip())
    return ids


def sample_of(path: str, sample_ids: list[str]) -> str | None:
    """Sample ID carried by a target path (bounded token match, longest ID first)."""
    name = Path(path).name
    for sid in sorted(sample_ids, key=len, reverse=True):
        if re.search(rf"(?:^|[_/]){re.escape(sid)}(?:[_.]|$)", name):
            return sid
    return None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Verify run completeness, excusing low-quality samples.")
    ap.add_argument("--results-dir", default="results", type=Path)
    ap.add_argument("--samplesheet", default="samplesheet.tsv", type=Path)
    ap.add_argument("--targets-file", required=True, type=Path,
                    help="File listing the run's required target paths, one per line.")
    ap.add_argument("--dest-name", default="low_quality_samples",
                    help="Subdirectory of results/ holding quarantined low-quality outputs.")
    args = ap.parse_args(argv)

    targets = read_lines(args.targets_file)
    if not targets:
        print("[verify] no required targets to check; treating run as complete.")
        return 0

    flagged = set(read_lines(args.results_dir / args.dest_name / "flagged_samples.txt"))
    sample_ids = read_sample_ids(args.samplesheet)

    missing = [t for t in targets if not Path(t).exists()]
    real_failures = [(t, sample_of(t, sample_ids)) for t in missing
                     if sample_of(t, sample_ids) not in flagged]
    excused = len(missing) - len(real_failures)

    if real_failures:
        print(f"[verify] RUN INCOMPLETE: {len(real_failures)} required output(s) missing for "
              f"sample(s) not flagged as low quality - marking the run FAILED:", file=sys.stderr)
        for target, sample in real_failures[:25]:
            print(f"  [{sample or 'unknown-sample'}] {target}", file=sys.stderr)
        if len(real_failures) > 25:
            print(f"  ... and {len(real_failures) - 25} more", file=sys.stderr)
        return 1

    print(f"[verify] run complete: all required outputs present for non-low-quality samples "
          f"({len(targets)} target(s) checked; {excused} missing output(s) excused for "
          f"{len(flagged)} quarantined low-quality sample(s)).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
