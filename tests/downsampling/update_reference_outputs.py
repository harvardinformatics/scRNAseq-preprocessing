import shutil
from pathlib import Path

from utils import TEST_SAMPLES, expected_output_for_sample, reference_output_for_sample, repo_root


def main():
    root = repo_root()
    copied = 0
    for sample in TEST_SAMPLES:
        source = root / expected_output_for_sample(sample)
        if not source.exists():
            raise SystemExit(f"missing workflow output: {source}")
        dest = root / reference_output_for_sample(sample)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, dest)
        copied += 1
    print(f"Copied {copied} reference output files")


if __name__ == "__main__":
    main()
