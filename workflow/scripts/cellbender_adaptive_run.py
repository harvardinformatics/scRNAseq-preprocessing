#!/usr/bin/env python
"""Adaptive CellBender run.

Run `cellbender remove-background`, read the automated assessment in the HTML report, and if it
recommends halving the learning rate, re-run at half the learning rate and keep whichever run's
learning curve is "normal" -- falling back to the initial run if the re-run is still not normal.
A per-sample status file records what happened, and BOTH runs' reports are kept for audit.

An unrecognized assessment (neither the "learning curve looks normal" nor the "re-run with half
the current learning rate" phrasing) is a HARD ERROR, so report-parsing edge cases surface
instead of being silently mishandled.

Invoked by workflow/rules/cellbender.smk inside the CellBender container. `assess_report` is a
pure function so it can be unit-tested on its own.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys

NORMAL_PHRASE = "this learning curve looks normal"
RERUN_PHRASE = "consider re-running with half the current learning rate"


def assess_report(html_text):
    """Classify a CellBender report's automated assessment.

    Returns 'normal', 'rerun', or 'unknown'. HTML tags are stripped and whitespace collapsed
    before matching, so minor formatting differences don't affect the result."""
    text = re.sub(r"<[^>]+>", " ", html_text)
    text = re.sub(r"\s+", " ", text).lower()
    if NORMAL_PHRASE in text:
        return "normal"
    if RERUN_PHRASE in text:
        return "rerun"
    return "unknown"


def assess_report_file(report_path):
    if not os.path.exists(report_path):
        raise FileNotFoundError(f"CellBender report not found: {report_path}")
    with open(report_path, encoding="utf-8", errors="replace") as handle:
        return assess_report(handle.read())


def cellbender_supports_seed():
    try:
        completed = subprocess.run(
            ["cellbender", "remove-background", "--help"],
            capture_output=True, text=True, check=False,
        )
    except FileNotFoundError:
        return False
    return "--seed" in (completed.stdout + completed.stderr)


def run_outputs(run_dir, base_name):
    """(base, filtered, report) paths CellBender writes for --output run_dir/base_name."""
    base = os.path.join(run_dir, base_name)
    stem = base[:-3] if base.endswith(".h5") else base
    return base, f"{stem}_filtered.h5", f"{stem}_report.html"


def run_cellbender(input_h5, output_h5, learning_rate, seed, seed_supported, run_dir):
    """Run cellbender remove-background in its own run_dir (own cwd -> own ckpt.tar.gz)."""
    os.makedirs(run_dir, exist_ok=True)
    cmd = [
        "cellbender", "remove-background", "--cuda",
        "--input", str(input_h5),
        "--output", str(output_h5),
        "--learning-rate", format(learning_rate, ".10g"),
    ]
    if seed_supported and seed is not None:
        cmd += ["--seed", str(seed)]
    print(f"[cellbender-adaptive] running: {' '.join(cmd)} (cwd={run_dir})", flush=True)
    subprocess.run(cmd, cwd=run_dir, check=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True)
    ap.add_argument("--output-base", required=True)
    ap.add_argument("--output-filtered", required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument("--status", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--learning-rate", type=float, required=True)
    ap.add_argument("--adaptive", required=True, help='"true" or "false"')
    ap.add_argument("--seed", default=None)
    ap.add_argument("--workdir", required=True)
    args = ap.parse_args(argv)

    adaptive = str(args.adaptive).strip().lower() in ("1", "true", "yes")
    input_h5 = os.path.abspath(args.input)
    output_base = os.path.abspath(args.output_base)
    output_filtered = os.path.abspath(args.output_filtered)
    output_report = os.path.abspath(args.report)
    output_status = os.path.abspath(args.status)
    base_name = os.path.basename(output_base)  # cellbender_{sample}.h5
    seed_supported = cellbender_supports_seed()

    workdir = os.path.abspath(args.workdir)
    initial_dir = os.path.join(workdir, "initial")
    rerun_dir = os.path.join(workdir, "rerun")
    initial_lr = args.learning_rate

    # --- initial run -------------------------------------------------------
    init_base, init_filtered, init_report = run_outputs(initial_dir, base_name)
    run_cellbender(input_h5, init_base, initial_lr, args.seed, seed_supported, initial_dir)
    initial_verdict = assess_report_file(init_report)
    print(f"[cellbender-adaptive] {args.sample}: initial run "
          f"(lr={format(initial_lr, '.10g')}) -> {initial_verdict}", flush=True)
    if initial_verdict == "unknown":
        sys.exit(
            f"[cellbender-adaptive] {args.sample}: could not classify the initial report "
            f"({init_report}) - found neither the 'learning curve looks normal' nor the "
            "'re-run with half the current learning rate' phrasing. Inspect the report and "
            "extend workflow/scripts/cellbender_adaptive_run.py to handle this case."
        )

    reran = False
    rerun_lr = None
    rerun_verdict = None
    chosen_dir = initial_dir

    if initial_verdict == "normal":
        outcome = "NORMAL_FIRST_TRY"
    elif not adaptive:
        # initial suggested a re-run, but adaptive mode is off -> keep the single run.
        outcome = "RERUN_SUGGESTED_BUT_ADAPTIVE_OFF"
    else:
        reran = True
        rerun_lr = initial_lr / 2.0
        rr_base, rr_filtered, rr_report = run_outputs(rerun_dir, base_name)
        run_cellbender(input_h5, rr_base, rerun_lr, args.seed, seed_supported, rerun_dir)
        rerun_verdict = assess_report_file(rr_report)
        print(f"[cellbender-adaptive] {args.sample}: re-run "
              f"(lr={format(rerun_lr, '.10g')}) -> {rerun_verdict}", flush=True)
        if rerun_verdict == "unknown":
            sys.exit(
                f"[cellbender-adaptive] {args.sample}: could not classify the re-run report "
                f"({rr_report}). Inspect the report and extend "
                "workflow/scripts/cellbender_adaptive_run.py to handle this case."
            )
        if rerun_verdict == "normal":
            chosen_dir = rerun_dir
            outcome = "RERUN_RESOLVED"
        else:
            chosen_dir = initial_dir
            outcome = "RERUN_DID_NOT_RESOLVE"

    kept_run = "rerun" if chosen_dir == rerun_dir else "initial"

    # --- copy the chosen run's outputs to the declared paths ---------------
    chosen_base, chosen_filtered, chosen_report = run_outputs(chosen_dir, base_name)
    os.makedirs(os.path.dirname(output_base), exist_ok=True)
    shutil.copy2(chosen_base, output_base)
    shutil.copy2(chosen_filtered, output_filtered)
    shutil.copy2(chosen_report, output_report)

    # --- keep BOTH runs' reports for audit ---------------------------------
    stem = output_report[:-len("_report.html")] if output_report.endswith("_report.html") \
        else os.path.splitext(output_report)[0]
    initial_archive = f"{stem}_report_initial.html"
    shutil.copy2(init_report, initial_archive)
    rerun_archive = "NA"
    if reran:
        rerun_archive = f"{stem}_report_rerun.html"
        shutil.copy2(rr_report, rerun_archive)

    # --- per-sample status -------------------------------------------------
    fields = [
        ("sample", args.sample),
        ("outcome", outcome),
        ("adaptive", str(adaptive).lower()),
        ("initial_learning_rate", format(initial_lr, ".10g")),
        ("initial_verdict", initial_verdict),
        ("reran", str(reran).lower()),
        ("rerun_learning_rate", format(rerun_lr, ".10g") if rerun_lr is not None else "NA"),
        ("rerun_verdict", rerun_verdict if rerun_verdict is not None else "NA"),
        ("kept_run", kept_run),
        ("report", output_report),
        ("report_initial", initial_archive),
        ("report_rerun", rerun_archive),
    ]
    with open(output_status, "w") as handle:
        for key, value in fields:
            handle.write(f"{key}\t{value}\n")
    print(f"[cellbender-adaptive] {args.sample}: outcome={outcome}, kept={kept_run}", flush=True)


if __name__ == "__main__":
    main()
