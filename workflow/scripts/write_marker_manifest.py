from pathlib import Path


def read_cluster_ids(cluster_ids_path):
    with open(cluster_ids_path) as handle:
        return [line.strip() for line in handle if line.strip()]


manifest_dir = Path(snakemake.output.manifest)
manifest_dir.mkdir(parents=True, exist_ok=True)

for existing_path in manifest_dir.iterdir():
    if existing_path.is_file():
        existing_path.unlink()

for cluster_id in read_cluster_ids(snakemake.input.cluster_ids):
    (manifest_dir / f"{cluster_id}.txt").write_text(f"{cluster_id}\\n")
