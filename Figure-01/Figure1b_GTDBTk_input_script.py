#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
from pathlib import Path


INPUT_CSV = Path("/scratch/p312334/project/10--HGT_isolates/data/__binstoget_isolates_summary_v3.csv")
OUTDIR = Path(os.environ.get("GTDBTK_OUTDIR", "/scratch/hb-fu/haoran/HGT_tree/gtdbtk_all_v3_rooted"))
GENOME_DIR = OUTDIR / "genomes"
META = OUTDIR / "meta_gtdbtk.txt"
ROOT_SOURCE = Path("/scratch/hb-fu/tkrempel/GTDB_tree/function25/genomes/root.fa")
ROOT_TAXONOMY = (
    "d__Bacteria;p__Minisyncoccota;c__Minisyncoccia;"
    "o__Minisyncoccales;f__Minisyncoccaceae;"
    "g__Minisyncoccus_A;s__Minisyncoccus archaeiphilus"
)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    GENOME_DIR.mkdir(parents=True, exist_ok=True)

    root_link = GENOME_DIR / "root.fa"
    if root_link.is_symlink() or root_link.exists():
        root_link.unlink()
    root_link.symlink_to(ROOT_SOURCE)

    seen = set()
    missing = []
    rows_written = 0

    with INPUT_CSV.open(newline="") as handle, META.open("w") as meta_handle:
        reader = csv.DictReader(handle)
        meta_handle.write(f"root\t{ROOT_TAXONOMY}\n")

        for row in reader:
            genome_id = row["Genome_file"].strip()
            fullpath = row["Fullpath"].strip()
            classification = row["classification"].strip()

            if not genome_id or genome_id in seen:
                continue
            seen.add(genome_id)

            if not fullpath or not classification:
                missing.append((genome_id, fullpath, classification, "missing_field"))
                continue

            src = Path(fullpath)
            if not src.exists():
                missing.append((genome_id, fullpath, classification, "missing_file"))
                continue

            link_path = GENOME_DIR / f"{genome_id}.fa"
            if link_path.is_symlink() or link_path.exists():
                link_path.unlink()
            link_path.symlink_to(src)

            meta_handle.write(f"{genome_id}\t{classification}\n")
            rows_written += 1

    missing_path = OUTDIR / "missing_or_skipped_genomes.tsv"
    with missing_path.open("w") as handle:
        handle.write("Genome_file\tFullpath\tclassification\treason\n")
        for genome_id, fullpath, classification, reason in missing:
            handle.write(f"{genome_id}\t{fullpath}\t{classification}\t{reason}\n")

    print(f"[INFO] Unique genome IDs written: {rows_written}")
    print(f"[INFO] Genome symlinks directory: {GENOME_DIR}")
    print(f"[INFO] Meta file: {META}")
    print(f"[INFO] Missing/skipped records: {len(missing)}")


if __name__ == "__main__":
    main()
