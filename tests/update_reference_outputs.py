#!/usr/bin/env python3
import shutil
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
MANIFEST = Path(__file__).with_name("test_reference_output_files.txt")
REFERENCE_ROOT = Path(__file__).with_name("reference_outputs")


def read_manifest(path):
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def main():
    outputs = read_manifest(MANIFEST)
    missing = [path for path in outputs if not (REPO_ROOT / path).exists()]
    if missing:
        raise SystemExit(
            "Cannot update reference outputs; missing current outputs:\n"
            + "\n".join(missing)
        )

    for rel_path in outputs:
        source = REPO_ROOT / rel_path
        target = REFERENCE_ROOT / rel_path
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)

    print(f"Copied {len(outputs)} reference output files to {REFERENCE_ROOT}")


if __name__ == "__main__":
    main()
