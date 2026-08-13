#!/usr/bin/env python3
"""Quarantine outputs of samples flagged as low-quality (too few cells to cluster).

The require_min_cells_for_pca() guard in workflow/scripts/silhouette_utils.R writes a
stable "[LOW_QUALITY_SAMPLE]" token into a job's log when QC/filtering leaves too few
cells to run PCA/clustering for a sample. This script scans the run logs for that token,
identifies the affected samples, and moves ALL of each flagged sample's output files from
results/<subdir>/ into results/low_quality_samples/<subdir>/, preserving the subdirectory
structure.

Job logs are deliberately left in results/logs/ (they explain why a sample was flagged,
and keeping them makes detection stable across re-runs). The move is per-sample: if a
sample is flagged in any branch/stage, every output carrying that sample name is moved.

The flagged sample IDs are also written to <results>/<dest>/flagged_samples.txt as an
ADVISORY record. This file does NOT drive the workflow. To actually skip a sample on
future runs, add it to `excluded_samples` in config/config.yaml (an explicit, visible
choice) - this script never edits the config, so a re-sequenced library reusing a sample
ID is never silently skipped.

Intended to run on completion of the Snakemake workflow (wired into the runner script).
Safe to run repeatedly; a no-op when nothing is flagged.
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

# Must match the token emitted by require_min_cells_for_pca() in silhouette_utils.R.
TOKEN = "[LOW_QUALITY_SAMPLE]"

# Top-level results subdirectories that are never scanned as outputs and never moved.
SKIP_DIRS = {"logs", "low_quality_samples", "low_quality_flags"}

# Log subdirectories skipped when scanning for the token, purely for speed. The guard runs
# only in the RunPCA scripts, so find_markers logs (the overwhelming majority - one per
# cluster per prefix) can never carry the token. Skipping them is safe: at worst, if these
# logs ever move, the scan just gets slower - it can never miss a real flag.
SKIP_LOG_DIRS = {"markers"}


def read_sample_ids(samplesheet: Path) -> list[str]:
    """Sample IDs from the first (tab-separated) column of the sample sheet, minus header."""
    ids: list[str] = []
    with samplesheet.open() as fh:
        next(fh, None)  # skip header row
        for line in fh:
            if not line.strip():
                continue
            ids.append(line.split("\t")[0].strip())
    return ids


def name_matches_sample(name: str, sid: str) -> bool:
    """True if a file/dir name carries the sample ID as a whole token.

    Sample names appear bounded by start/underscore/slash on the left and by
    underscore/dot/end on the right, e.g. "..._cteleta.rds", "cteleta_..._matrix",
    "cellbender_cteleta.h5". Bounded matching prevents one sample ID from matching
    inside another (e.g. sc108 vs sc109, obOLsample1 vs obOLsample2).
    """
    return re.search(rf"(?:^|[_/]){re.escape(sid)}(?:[_.]|$)", name) is not None


def _token_log_names(logs_dir: Path) -> list[str]:
    """Basenames of log files containing the token. Uses grep (fast over a large, networked
    log tree - this workflow can produce tens of thousands of logs) and falls back to a
    pure-Python scan if grep is unavailable."""
    exclude_args = [f"--exclude-dir={d}" for d in SKIP_LOG_DIRS]
    try:
        result = subprocess.run(
            ["grep", "-rlF", *exclude_args, TOKEN, str(logs_dir)],
            capture_output=True, text=True, check=False,
        )
        # rc 0 = matches, 1 = no matches; anything else -> fall back.
        if result.returncode in (0, 1):
            return [Path(line).name for line in result.stdout.splitlines() if line]
    except (FileNotFoundError, OSError):
        pass
    names: list[str] = []
    for log_path in logs_dir.rglob("*.log"):
        if any(part in SKIP_LOG_DIRS for part in log_path.relative_to(logs_dir).parts):
            continue
        try:
            if TOKEN in log_path.read_text(errors="replace"):
                names.append(log_path.name)
        except OSError:
            continue
    return names


def find_flagged_samples(results_dir: Path, sample_ids: list[str]) -> set[str]:
    logs_dir = results_dir / "logs"
    flagged: set[str] = set()
    if not logs_dir.is_dir():
        return flagged
    # Longest IDs first so the most specific sample name wins the attribution.
    ordered = sorted(sample_ids, key=len, reverse=True)
    for name in _token_log_names(logs_dir):
        for sid in ordered:
            if name_matches_sample(name, sid):
                flagged.add(sid)
                break
    return flagged


def quarantine(results_dir: Path, flagged: set[str], dest_name: str,
               dry_run: bool) -> list[tuple[Path, Path]]:
    dest_root = results_dir / dest_name
    moves: list[tuple[Path, Path]] = []
    for sub in sorted(p for p in results_dir.iterdir() if p.is_dir()):
        if sub.name in SKIP_DIRS:
            continue
        for item in sorted(sub.iterdir()):
            if not any(name_matches_sample(item.name, sid) for sid in flagged):
                continue
            dest_dir = dest_root / sub.name
            dest = dest_dir / item.name
            moves.append((item, dest))
            if dry_run:
                continue
            dest_dir.mkdir(parents=True, exist_ok=True)
            if dest.is_dir():
                shutil.rmtree(dest)
            elif dest.exists():
                dest.unlink()
            shutil.move(str(item), str(dest))
    return moves


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Quarantine low-quality-sample outputs.")
    ap.add_argument("--results-dir", default="results", type=Path,
                    help="Workflow results directory (default: results).")
    ap.add_argument("--samplesheet", default="samplesheet.tsv", type=Path,
                    help="Sample sheet whose first column lists sample IDs (default: samplesheet.tsv).")
    ap.add_argument("--dest-name", default="low_quality_samples",
                    help="Subdirectory of results/ to move flagged outputs into.")
    ap.add_argument("--dry-run", action="store_true",
                    help="Report what would move without moving anything.")
    args = ap.parse_args(argv)

    if not args.results_dir.is_dir():
        print(f"[quarantine] results dir not found: {args.results_dir}; nothing to do.")
        return 0
    if not args.samplesheet.is_file():
        print(f"[quarantine] sample sheet not found: {args.samplesheet}; skipping quarantine.",
              file=sys.stderr)
        return 0

    sample_ids = read_sample_ids(args.samplesheet)
    flagged = find_flagged_samples(args.results_dir, sample_ids)
    if not flagged:
        print("[quarantine] no low-quality samples flagged; nothing to move.")
        return 0

    flagged_sorted = sorted(flagged)
    print(f"[quarantine] flagged low-quality sample(s): {', '.join(flagged_sorted)}")

    # Write the flagged manifest BEFORE moving anything. It records what was detected and is
    # also read by verify_run_complete.py to excuse these samples from the run's pass/fail
    # decision; writing it first keeps that decision correct even if a later move errors.
    dest_root = args.results_dir / args.dest_name
    manifest = dest_root / "flagged_samples.txt"
    if not args.dry_run:
        dest_root.mkdir(parents=True, exist_ok=True)
        manifest.write_text("\n".join(flagged_sorted) + "\n")

    moves = quarantine(args.results_dir, flagged, args.dest_name, args.dry_run)
    verb = "would move" if args.dry_run else "moved"
    for src, dest in moves:
        print(f"[quarantine] {verb}: {src} -> {dest}")
    print(f"[quarantine] {verb} {len(moves)} item(s) for {len(flagged)} flagged sample(s) "
          f"into {args.results_dir / args.dest_name}/")

    # The manifest is advisory for the DAG: it does NOT exclude samples. To skip these on
    # future runs, add them to excluded_samples in config/config.yaml (an explicit choice).
    print(f"[quarantine] advisory list {'would be ' if args.dry_run else ''}written to {manifest}")
    print("[quarantine] to skip these samples on future runs, add them to "
          "'excluded_samples' in config/config.yaml")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
