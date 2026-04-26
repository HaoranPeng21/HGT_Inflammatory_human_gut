#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import hashlib
import re
import pandas as pd


BASE = Path("/scratch/p312334/project/10--HGT_isolates/2--Phenotype_association/comparision_plot/Streptococcus parasanguinis D | Bifidobacterium longum")
RAW = BASE / "233_3-338_1_raw_HGT_events.csv"
PANGENOME = Path("/scratch/p312334/project/10--HGT_isolates/pangenome")
OUT = Path("/scratch/p312334/project/10--HGT_isolates/2--Phenotype_association/comparision_plot/233_3-338_1_plasmid_exactdedup_lovis4u_v2")
GFF_OUT = OUT / "gff_input"


def sanitize_name(text: str) -> str:
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text)
    return re.sub(r"_+", "_", text).strip("_")


def find_gff(genome: str) -> Path:
    matches = sorted(PANGENOME.glob(f"**/fixed_input_files/{genome}.unicycler.gff"))
    if not matches:
        raise FileNotFoundError(f"GFF not found for genome {genome}")
    return matches[0]


def read_gff(path: Path) -> list[str]:
    return path.read_text(encoding="utf-8", errors="replace").splitlines()


def extract_single_contig(lines: list[str], contig: str) -> tuple[list[str], str]:
    out_lines: list[str] = []
    in_fasta = False
    keep_fasta = False
    seq_chunks: list[str] = []

    for line in lines:
        if line.startswith("##FASTA"):
            in_fasta = True
            out_lines.append(line)
            continue
        if not in_fasta:
            if line.startswith("#"):
                out_lines.append(line)
                continue
            if not line.strip():
                continue
            seqid = line.split("\t", 1)[0]
            if seqid == contig:
                out_lines.append(line)
        else:
            if line.startswith(">"):
                current_id = line[1:].split()[0]
                keep_fasta = current_id == contig
                if keep_fasta:
                    out_lines.append(line)
            elif keep_fasta:
                out_lines.append(line)
                seq_chunks.append(line.strip())

    seq = "".join(seq_chunks).upper()
    if not seq:
        raise ValueError(f"No FASTA sequence found for contig {contig}")
    return out_lines, seq


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    GFF_OUT.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(RAW)
    candidates: list[dict[str, str]] = []

    for _, row in df.iterrows():
        if int(row["plasmid1"]) == 1:
            candidates.append(
                {
                    "side": "query",
                    "genome": str(row["MAG_1"]),
                    "contig_id": str(row["query_contig_id"]),
                }
            )
        if int(row["plamids2"]) == 1:
            candidates.append(
                {
                    "side": "subject",
                    "genome": str(row["MAG_2"]),
                    "contig_id": str(row["subject_contig_id"]),
                }
            )

    cand_df = pd.DataFrame(candidates).drop_duplicates().sort_values(["side", "genome", "contig_id"])

    raw_rows = []
    for _, row in cand_df.iterrows():
        genome = row["genome"]
        contig = row["contig_id"]
        side = row["side"]
        gff = find_gff(genome)
        lines = read_gff(gff)
        subset_lines, seq = extract_single_contig(lines, contig)
        seq_hash = hashlib.sha1(seq.encode()).hexdigest()
        raw_rows.append(
            {
                "side": side,
                "genome": genome,
                "contig_id": contig,
                "source_gff": str(gff),
                "length": len(seq),
                "seq_hash": seq_hash,
                "subset_lines": subset_lines,
                "sequence": seq,
            }
        )

    full_df = pd.DataFrame(raw_rows)
    full_df["length_ratio_to_first"] = full_df.groupby("seq_hash")["length"].transform(lambda s: s / s.iloc[0])

    dedup_groups = []
    rep_rows = []
    groups = []
    for seq_hash, sub in full_df.groupby("seq_hash", sort=False):
        groups.append((seq_hash, sub.copy()))

    groups.sort(key=lambda x: (-int(x[1]["length"].iloc[0]), x[1]["genome"].iloc[0], x[1]["contig_id"].iloc[0]))

    for i, (seq_hash, sub) in enumerate(groups, start=1):
        sub = sub.sort_values(["side", "genome", "contig_id"]).copy()
        # User rule: keep one when sequences are identical and lengths are within 90%.
        # Exact sequence identity already implies equal length, so this condition is satisfied.
        rep = sub.iloc[0]
        label = f"{i:02d}_Module_{chr(64+i)}_{rep['length']}bp_n{len(sub)}"
        out_gff = GFF_OUT / f"{label}.gff"
        out_gff.write_text("\n".join(rep["subset_lines"]) + "\n", encoding="utf-8")

        members = ";".join(f"{r.side}:{r.genome}:c{r.contig_id}" for r in sub.itertuples())
        rep_rows.append(
            {
                "rep_label": label,
                "rep_side": rep["side"],
                "rep_genome": rep["genome"],
                "rep_contig_id": rep["contig_id"],
                "length": rep["length"],
                "n_members": len(sub),
                "members": members,
                "output_gff": str(out_gff),
                "seq_hash": seq_hash,
            }
        )
        for r in sub.itertuples():
            dedup_groups.append(
                {
                    "rep_label": label,
                    "side": r.side,
                    "genome": r.genome,
                    "contig_id": r.contig_id,
                    "length": r.length,
                    "seq_hash": seq_hash,
                }
            )

    pd.DataFrame(rep_rows).to_csv(OUT / "representative_manifest.csv", index=False)
    pd.DataFrame(dedup_groups).to_csv(OUT / "dedup_membership.csv", index=False)

    summary = [
        f"candidate_plasmid_contigs\t{len(full_df)}",
        f"representative_contigs\t{len(rep_rows)}",
        f"query_plasmid_candidates\t{(full_df['side'] == 'query').sum()}",
        f"subject_plasmid_candidates\t{(full_df['side'] == 'subject').sum()}",
    ]
    (OUT / "summary.tsv").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
