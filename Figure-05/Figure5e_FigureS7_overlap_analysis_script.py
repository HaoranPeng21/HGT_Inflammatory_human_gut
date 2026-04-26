#!/usr/bin/env python3
import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

BASE = Path("/scratch/p312334/project")
PYSEER_BASE = BASE / "10--HGT_isolates" / "3.2--Pangenome"
PANGENOME_BASE = BASE / "10--HGT_isolates" / "pangenome"
HGT_FILE = BASE / "10--HGT_isolates" / "data" / "HGT_table_2000bp_isolates_within_plasmid.csv"
TOP_N = 50
FDR_THRESH = 0.05
ATTR_RE = re.compile(r"([^=;]+)=([^;]*)")
SKIP_DIRS = {"jobs", "log", "venv", "__pycache__"}


def norm_genome(name: str) -> str:
    text = str(name)
    return text[:-10] if text.endswith(".unicycler") else text


def norm_contig(value) -> str:
    text = str(value).strip()
    if text.endswith(".0"):
        text = text[:-2]
    return text


def sort_interval(start, end):
    a = int(float(start))
    b = int(float(end))
    return (a, b) if a <= b else (b, a)


def overlap_len(a_start, a_end, b_start, b_end):
    left = max(a_start, b_start)
    right = min(a_end, b_end)
    return max(0, right - left + 1)


def bh_fdr(pvals):
    n = len(pvals)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvals[i])
    qvals = [1.0] * n
    prev = 1.0
    for pos, idx in enumerate(reversed(order), start=1):
        rank = n - pos + 1
        q = min(prev, pvals[idx] * n / rank)
        qvals[idx] = q
        prev = q
    return qvals


def find_species_dirs(base: Path):
    species_dirs = []
    for path in sorted(base.iterdir()):
        if not path.is_dir() or path.name in SKIP_DIRS or path.name.startswith("_"):
            continue
        if (path / "IBD_COGs.txt").exists():
            species_dirs.append(path)
    return species_dirs


def read_pyseer(pyseer_file: Path):
    rows = []
    with pyseer_file.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                p = float(row["lrt-pvalue"])
                beta = float(row["beta"])
                af = float(row["af"])
            except Exception:
                continue
            if p <= 0:
                continue
            row["_p"] = p
            row["_beta"] = beta
            row["_af"] = af
            row["_logp"] = -math.log10(p)
            rows.append(row)
    signal_rows = [r for r in rows if r.get("notes", "").strip() == ""]
    qvals = bh_fdr([r["_p"] for r in signal_rows])
    for row, q in zip(signal_rows, qvals):
        row["_fdr"] = q
    sig = [r for r in signal_rows if r["_fdr"] < FDR_THRESH]
    sig.sort(key=lambda r: (r["_fdr"], r["_p"], -abs(r["_beta"])))
    top = sig[:TOP_N]
    return rows, signal_rows, sig, top


def write_sig_tables(pyseer_dir: Path, sig, top):
    sig_file = pyseer_dir / "IBD_COGs_FDR_lt_0.05.tsv"
    top_file = pyseer_dir / "IBD_COGs_FDR_lt_0.05_top50.tsv"
    fields = ["variant", "af", "lrt-pvalue", "beta", "notes", "fdr", "log10p"]
    with sig_file.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in sig:
            w.writerow({
                "variant": r["variant"],
                "af": r["af"],
                "lrt-pvalue": r["lrt-pvalue"],
                "beta": r["beta"],
                "notes": r.get("notes", ""),
                "fdr": f"{r['_fdr']:.6g}",
                "log10p": f"{r['_logp']:.6f}",
            })
    with top_file.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in top:
            w.writerow({
                "variant": r["variant"],
                "af": r["af"],
                "lrt-pvalue": r["lrt-pvalue"],
                "beta": r["beta"],
                "notes": r.get("notes", ""),
                "fdr": f"{r['_fdr']:.6g}",
                "log10p": f"{r['_logp']:.6f}",
            })
    return sig_file, top_file


def load_hgt_index(target_mags):
    hgt_index = defaultdict(list)
    matched_mags = set()
    with HGT_FILE.open(newline="") as f:
        reader = csv.DictReader(f)
        for row_num, row in enumerate(reader, start=2):
            q_start, q_end = sort_interval(row["q_start"], row["q_end"])
            s_start, s_end = sort_interval(row["s_start"], row["s_end"])
            mag1 = norm_genome(row["MAG_1"])
            mag2 = norm_genome(row["MAG_2"])
            if mag1 not in target_mags and mag2 not in target_mags:
                continue
            contig1 = norm_contig(row["query_contig_id"])
            contig2 = norm_contig(row["subject_contig_id"])
            common = {
                "row_num": row_num,
                "MAG_1": row["MAG_1"],
                "MAG_2": row["MAG_2"],
                "Sequence_query": row["Sequence_query"],
                "Sequence_subject": row["Sequence_subject"],
                "species_pair": row["species_pair"],
                "individual_1": row["individual_1"],
                "individual_2": row["individual_2"],
                "Person": row["Person"],
                "within_cohort": row["within_cohort"],
                "Identity": row["Identity"],
                "Length": row["Length"],
                "plasmid1": row.get("plasmid1", ""),
                "plamids2": row.get("plamids2", ""),
                "plasmid": row.get("plasmid", ""),
            }
            if mag1 in target_mags:
                matched_mags.add(mag1)
                hgt_index[(mag1, contig1)].append({
                    **common,
                    "hgt_side": "MAG_1",
                    "hgt_mag": mag1,
                    "hgt_contig": contig1,
                    "hgt_start": q_start,
                    "hgt_end": q_end,
                })
            if mag2 in target_mags:
                matched_mags.add(mag2)
                hgt_index[(mag2, contig2)].append({
                    **common,
                    "hgt_side": "MAG_2",
                    "hgt_mag": mag2,
                    "hgt_contig": contig2,
                    "hgt_start": s_start,
                    "hgt_end": s_end,
                })
    return hgt_index, matched_mags


def load_roary_rows(roary_file: Path, selected, allowed_mags):
    roary_rows = {}
    with roary_file.open(newline="") as f:
        reader = csv.DictReader(f)
        genome_cols = [col.strip('"') for col in reader.fieldnames[14:]]
        for row in reader:
            gene = row["Gene"].strip('"')
            if gene not in selected:
                continue
            genomes = {}
            for col in genome_cols:
                raw_value = row.get(col, "")
                if raw_value and norm_genome(col) in allowed_mags:
                    genomes[col] = raw_value.strip('"')
            if genomes:
                roary_rows[gene] = {
                    "annotation": row["Annotation"].strip('"'),
                    "genomes": genomes,
                }
    return roary_rows


def parse_coords(gff_dir: Path, roary_rows):
    needed = defaultdict(set)
    for info in roary_rows.values():
        for genome, gene_ids in info["genomes"].items():
            for gene_id in gene_ids.split(";"):
                gene_id = gene_id.strip()
                if gene_id:
                    needed[genome].add(gene_id)

    coords = {}
    for genome, ids in needed.items():
        gff = gff_dir / f"{genome}.gff"
        if not gff.exists():
            gff = gff_dir / f"{genome}.gff3"
        if not gff.exists():
            continue
        with gff.open() as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9 or parts[2] != "CDS":
                    continue
                attrs = dict(ATTR_RE.findall(parts[8]))
                gene_id = attrs.get("ID") or attrs.get("locus_tag")
                if gene_id not in ids:
                    continue
                coords[(genome, gene_id)] = {
                    "contig": norm_contig(parts[0]),
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "strand": parts[6],
                    "gene_name": attrs.get("gene", ""),
                    "product": attrs.get("product", ""),
                }
    return coords


def write_overlap(pyseer_dir: Path, roary_rows, coords, hgt_index, top):
    detail_file = pyseer_dir / "top50_FDR_lt_0.05_within_HGT_overlaps.tsv"
    summary_file = pyseer_dir / "top50_FDR_lt_0.05_within_HGT_overlap_summary.tsv"
    stats = {r["variant"]: r for r in top}
    detail_rows = []
    by_group = defaultdict(lambda: {"gene_instances": set(), "hgt_rows": set(), "max_overlap_bp": 0})

    for group, info in roary_rows.items():
        for genome, gene_ids in info["genomes"].items():
            ngenome = norm_genome(genome)
            for gene_id in gene_ids.split(";"):
                gene_id = gene_id.strip()
                coord = coords.get((genome, gene_id))
                if not coord:
                    continue
                hits = hgt_index.get((ngenome, coord["contig"]), [])
                if not hits:
                    continue
                gene_start, gene_end = sort_interval(coord["start"], coord["end"])
                gene_len = gene_end - gene_start + 1
                for hit in hits:
                    ov = overlap_len(gene_start, gene_end, hit["hgt_start"], hit["hgt_end"])
                    if ov == 0:
                        continue
                    hgt_len = hit["hgt_end"] - hit["hgt_start"] + 1
                    detail_rows.append({
                        "group": group,
                        "group_annotation": info["annotation"],
                        "lrt_pvalue": stats[group]["lrt-pvalue"],
                        "fdr": f"{stats[group]['_fdr']:.6g}",
                        "log10p": f"{stats[group]['_logp']:.6f}",
                        "beta": stats[group]["beta"],
                        "af": stats[group]["af"],
                        "genome": genome,
                        "gene_id": gene_id,
                        "contig": coord["contig"],
                        "start": gene_start,
                        "end": gene_end,
                        "strand": coord["strand"],
                        "gene_name": coord["gene_name"],
                        "product": coord["product"],
                        "row_num": hit["row_num"],
                        "hgt_side": hit["hgt_side"],
                        "hgt_mag": hit["hgt_mag"],
                        "hgt_contig": hit["hgt_contig"],
                        "hgt_start": hit["hgt_start"],
                        "hgt_end": hit["hgt_end"],
                        "overlap_bp": ov,
                        "gene_overlap_fraction": f"{ov / gene_len:.6f}",
                        "hgt_overlap_fraction": f"{ov / hgt_len:.6f}",
                        "MAG_1": hit["MAG_1"],
                        "MAG_2": hit["MAG_2"],
                        "Sequence_query": hit["Sequence_query"],
                        "Sequence_subject": hit["Sequence_subject"],
                        "species_pair": hit["species_pair"],
                        "individual_1": hit["individual_1"],
                        "individual_2": hit["individual_2"],
                        "Person": hit["Person"],
                        "within_cohort": hit["within_cohort"],
                        "Identity": hit["Identity"],
                        "Length": hit["Length"],
                        "plasmid1": hit["plasmid1"],
                        "plamids2": hit["plamids2"],
                        "plasmid": hit["plasmid"],
                    })
                    by_group[group]["gene_instances"].add((genome, gene_id))
                    by_group[group]["hgt_rows"].add(hit["row_num"])
                    by_group[group]["max_overlap_bp"] = max(by_group[group]["max_overlap_bp"], ov)

    if detail_rows:
        with detail_file.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(detail_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(detail_rows)
    else:
        detail_file.write_text("")

    summary_rows = []
    for group in [r["variant"] for r in top]:
        if group not in by_group:
            continue
        s = stats[group]
        info = by_group[group]
        summary_rows.append({
            "group": group,
            "group_annotation": roary_rows.get(group, {}).get("annotation", ""),
            "lrt_pvalue": s["lrt-pvalue"],
            "fdr": f"{s['_fdr']:.6g}",
            "log10p": f"{s['_logp']:.6f}",
            "beta": s["beta"],
            "af": s["af"],
            "gene_instances": len(info["gene_instances"]),
            "hgt_rows": len(info["hgt_rows"]),
            "max_overlap_bp": info["max_overlap_bp"],
        })
    if summary_rows:
        with summary_file.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(summary_rows)
    else:
        summary_file.write_text("")
    return detail_file, summary_file, summary_rows


def make_plot(species: str, pyseer_dir: Path, all_rows, summary_rows):
    plot_file = pyseer_dir / "IBD_COGs_top50_FDR_HGT_highlight.png"
    highlight = {r["group"] for r in summary_rows}
    hrows = [r for r in all_rows if r["variant"] in highlight]
    plt.figure(figsize=(10, 7))
    sc = plt.scatter(
        [r["_beta"] for r in all_rows],
        [r["_logp"] for r in all_rows],
        c=[r["_af"] for r in all_rows],
        cmap="viridis",
        s=14,
        alpha=0.65,
        edgecolors="none",
    )
    if hrows:
        plt.scatter(
            [r["_beta"] for r in hrows],
            [r["_logp"] for r in hrows],
            facecolors="none",
            edgecolors="#d62728",
            s=70,
            linewidths=1.0,
        )
        for r in hrows:
            plt.text(r["_beta"], r["_logp"], r["variant"], fontsize=7, ha="left", va="bottom")
    cbar = plt.colorbar(sc)
    cbar.set_label("Allele frequency")
    plt.xlabel("Effect size (beta)")
    plt.ylabel("-log10(p-value)")
    plt.title(f"{species} IBD pyseer\nTop50 FDR<0.05 with within-HGT hits highlighted")
    plt.tight_layout()
    plt.savefig(plot_file, dpi=300)
    plt.close()
    return plot_file


def analyze_species(species: str):
    pyseer_dir = PYSEER_BASE / species
    pangenome_dir = PANGENOME_BASE / species
    pyseer_file = pyseer_dir / "IBD_COGs.txt"
    roary_file = pangenome_dir / "gene_presence_absence.csv"
    gff_dir = pangenome_dir / "fixed_input_files"
    if not pyseer_file.exists():
        raise FileNotFoundError(pyseer_file)
    if not roary_file.exists():
        raise FileNotFoundError(roary_file)
    if not gff_dir.exists():
        raise FileNotFoundError(gff_dir)

    all_rows, signal_rows, sig, top = read_pyseer(pyseer_file)
    sig_file, top_file = write_sig_tables(pyseer_dir, sig, top)
    with roary_file.open(newline="") as f:
        reader = csv.DictReader(f)
        genome_cols = [col.strip('"') for col in reader.fieldnames[14:]]
        target_mags = {norm_genome(c) for c in genome_cols}
    hgt_index, matched_mags = load_hgt_index(target_mags)
    roary_rows = load_roary_rows(roary_file, {r["variant"] for r in top}, matched_mags)
    coords = parse_coords(gff_dir, roary_rows)
    detail_file, summary_file, summary_rows = write_overlap(pyseer_dir, roary_rows, coords, hgt_index, top)
    plot_file = make_plot(species, pyseer_dir, all_rows, summary_rows)
    return {
        "species": species,
        "signal_rows": len(signal_rows),
        "sig_fdr_lt_0_05": len(sig),
        "top_selected": len(top),
        "matched_hgt_mags": len(matched_mags),
        "hgt_highlighted": len(summary_rows),
        "sig_file": str(sig_file),
        "top_file": str(top_file),
        "detail_file": str(detail_file),
        "summary_file": str(summary_file),
        "plot_file": str(plot_file),
    }


def main(argv):
    if len(argv) > 1:
        species_list = argv[1:]
    else:
        species_list = [p.name for p in find_species_dirs(PYSEER_BASE)]

    results = []
    for species in species_list:
        results.append(analyze_species(species))

    summary_path = PYSEER_BASE / "top50_FDR_lt_0.05_within_HGT_species_summary.tsv"
    fields = [
        "species", "signal_rows", "sig_fdr_lt_0_05", "top_selected",
        "matched_hgt_mags", "hgt_highlighted", "sig_file", "top_file",
        "detail_file", "summary_file", "plot_file",
    ]
    with summary_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(results)
    print(summary_path)
    for row in results:
        print(f"{row['species']}\tFDR<0.05={row['sig_fdr_lt_0_05']}\ttop50={row['top_selected']}\tHGT_hits={row['hgt_highlighted']}")


if __name__ == "__main__":
    main(sys.argv)
