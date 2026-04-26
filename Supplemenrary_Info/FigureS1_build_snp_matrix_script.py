#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np


BASE = Path("/scratch/p312334/project/10--HGT_isolates/1--General_infomation/conserve_cutoff")
UCG_DIR = BASE / "secondary_cluster_best.tmp"
ALIGN_DIR = BASE / "secondary_cluster_best.res" / "secondary_cluster_best"
OUTDIR = BASE / "secondary_cluster_best.res" / "ribosomal42_snp"

RIBO = [
    "rplA", "rplB", "rplC", "rplD", "rplE", "rplF", "rplI", "rplJ", "rplK", "rplL",
    "rplM", "rplN", "rplO", "rplP", "rplQ", "rplR", "rplS", "rplT", "rplU", "rplV",
    "rplW", "rplX", "rpmA", "rpmC", "rpmI", "rpsB", "rpsC", "rpsD", "rpsE", "rpsF",
    "rpsG", "rpsH", "rpsI", "rpsJ", "rpsL", "rpsM", "rpsO", "rpsP", "rpsQ", "rpsR",
    "rpsS", "rpsT",
]


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    current = None
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:]
                seqs[current] = []
            else:
                if current is None:
                    raise ValueError(f"Sequence data before header in {path}")
                seqs[current].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def load_uid_to_label() -> dict[str, str]:
    mapping = {}
    for path in sorted(UCG_DIR.glob("*.ucg")):
        with path.open() as handle:
            obj = json.load(handle)
        info = obj.get("genome_info", {})
        uid = info.get("uid")
        label = info.get("label")
        if uid is None or not label:
            continue
        mapping[f"zZ{uid}zZ"] = label
    return mapping


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    uid_to_label = load_uid_to_label()
    all_headers = sorted(uid_to_label)
    concatenated: dict[str, list[str]] = {h: [] for h in all_headers}
    marker_lengths: list[tuple[str, int]] = []

    for marker in RIBO:
        path = ALIGN_DIR / f"aligned_{marker}_codon.zZ.fasta"
        seqs = read_fasta(path)
        length_set = {len(v) for v in seqs.values()}
        if len(length_set) != 1:
            raise SystemExit(f"Unequal alignment lengths in {path}")
        aln_len = next(iter(length_set))
        marker_lengths.append((marker, aln_len))
        gap_seq = "-" * aln_len

        for header in all_headers:
            concatenated[header].append(seqs.get(header, gap_seq).upper())

    labels = [uid_to_label.get(h, h) for h in all_headers]
    seq_strings = ["".join(concatenated[h]) for h in all_headers]

    concat_fasta = OUTDIR / "ribosomal42_concatenated_codon.fasta"
    with concat_fasta.open("w") as handle:
        for label, seq in zip(labels, seq_strings):
            handle.write(f">{label}\n{seq}\n")

    lengths_tsv = OUTDIR / "ribosomal42_marker_lengths.tsv"
    with lengths_tsv.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["marker", "aligned_length_bp"])
        writer.writerows(marker_lengths)

    arr = np.frombuffer("".join(seq_strings).encode("ascii"), dtype=np.uint8).reshape(len(seq_strings), -1)
    valid = np.isin(arr, np.frombuffer(b"ACGT", dtype=np.uint8))

    n = arr.shape[0]
    snp_mat = np.zeros((n, n), dtype=np.int32)
    comparable_mat = np.zeros((n, n), dtype=np.int32)

    for i in range(n):
        valid_pair = valid[i] & valid[i:]
        diff = (arr[i] != arr[i:]) & valid_pair
        snp_counts = diff.sum(axis=1, dtype=np.int32)
        comparable = valid_pair.sum(axis=1, dtype=np.int32)
        snp_mat[i, i:] = snp_counts
        snp_mat[i:, i] = snp_counts
        comparable_mat[i, i:] = comparable
        comparable_mat[i:, i] = comparable

    matrix_tsv = OUTDIR / "ribosomal42_pairwise_snp_matrix.tsv"
    with matrix_tsv.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["genome"] + labels)
        for label, row in zip(labels, snp_mat):
            writer.writerow([label] + row.tolist())

    comparable_tsv = OUTDIR / "ribosomal42_pairwise_comparable_sites_matrix.tsv"
    with comparable_tsv.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["genome"] + labels)
        for label, row in zip(labels, comparable_mat):
            writer.writerow([label] + row.tolist())

    long_tsv = OUTDIR / "ribosomal42_pairwise_snp_long.tsv"
    with long_tsv.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["genome_1", "genome_2", "snp_count", "comparable_sites", "snp_rate", "snp_per_2000bp"])
        for i in range(n):
            for j in range(i + 1, n):
                comparable = int(comparable_mat[i, j])
                snps = int(snp_mat[i, j])
                rate = (snps / comparable) if comparable else ""
                per_2000 = (snps / comparable * 2000) if comparable else ""
                writer.writerow([labels[i], labels[j], snps, comparable, rate, per_2000])

    print(f"[INFO] Genomes: {n}")
    print(f"[INFO] Markers: {len(RIBO)}")
    print(f"[INFO] Concatenated length: {arr.shape[1]} bp")
    print(f"[INFO] Wrote: {concat_fasta}")
    print(f"[INFO] Wrote: {matrix_tsv}")
    print(f"[INFO] Wrote: {comparable_tsv}")
    print(f"[INFO] Wrote: {long_tsv}")


if __name__ == "__main__":
    main()
