#!/usr/bin/env python3
import csv
from collections import defaultdict
from pathlib import Path


BASE = Path("/scratch/p312334/project/10--HGT_isolates")
PYSEER_BASE = BASE / "3.2--Pangenome"
PANGENOME_BASE = BASE / "pangenome"
PYSEER_INPUT_BASE = BASE / "pangenome" / "pyseer_skani9999"
HGT_FILE = BASE / "data" / "HGT_table_2000bp_isolates_within_plasmid.csv"
FDR_THRESH = 0.05
MIN_SKANI9999_KEPT = 100


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


def gff_attr_id(attrs: str):
    if attrs.startswith("ID="):
        return attrs[3:].split(";", 1)[0]
    marker = ";ID="
    pos = attrs.find(marker)
    if pos >= 0:
        return attrs[pos + len(marker):].split(";", 1)[0]
    if attrs.startswith("locus_tag="):
        return attrs[10:].split(";", 1)[0]
    marker = ";locus_tag="
    pos = attrs.find(marker)
    if pos >= 0:
        return attrs[pos + len(marker):].split(";", 1)[0]
    return None


def species_result_dir(species: str) -> Path:
    direct = PYSEER_BASE / species
    if (direct / "IBD_COGs_FDR_lt_0.05.tsv").exists():
        return direct
    nested = PYSEER_BASE / "None" / species
    if (nested / "IBD_COGs_FDR_lt_0.05.tsv").exists():
        return nested
    raise FileNotFoundError(f"missing 3.2 result directory for {species}")


def read_input_summary(species: str):
    path = PYSEER_INPUT_BASE / species / "input_summary.tsv"
    data = {}
    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) >= 2:
                data[row[0]] = row[1]
    return data


def selected_species():
    species = []
    for path in sorted(PYSEER_INPUT_BASE.iterdir()):
        if not path.is_dir():
            continue
        summary = path / "input_summary.tsv"
        if not summary.exists():
            continue
        data = read_input_summary(path.name)
        kept = int(data.get("kept_samples", "0"))
        if kept <= MIN_SKANI9999_KEPT:
            continue
        try:
            species_result_dir(path.name)
        except FileNotFoundError:
            continue
        species.append((path.name, kept))
    return species


def read_sig_genes(species: str):
    sig_file = species_result_dir(species) / "IBD_COGs_FDR_lt_0.05.tsv"
    genes = []
    with sig_file.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                fdr = float(row["fdr"])
                beta = float(row["beta"])
            except (TypeError, ValueError):
                continue
            if fdr >= FDR_THRESH:
                continue
            row["_fdr"] = fdr
            row["_beta"] = beta
            genes.append(row)
    return genes


def read_genome_columns(roary_file: Path):
    with roary_file.open(newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
    return [col.strip('"') for col in header[14:]]


def load_hgt_index(target_mags):
    hgt_index = defaultdict(list)
    target_mags = set(target_mags)
    with HGT_FILE.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row_num, row in enumerate(reader, start=2):
            mag1 = norm_genome(row["MAG_1"])
            mag2 = norm_genome(row["MAG_2"])
            if mag1 not in target_mags and mag2 not in target_mags:
                continue
            q_start, q_end = sort_interval(row["q_start"], row["q_end"])
            s_start, s_end = sort_interval(row["s_start"], row["s_end"])
            if mag1 in target_mags:
                hgt_index[(mag1, norm_contig(row["query_contig_id"]))].append(
                    {"row_num": row_num, "start": q_start, "end": q_end}
                )
            if mag2 in target_mags:
                hgt_index[(mag2, norm_contig(row["subject_contig_id"]))].append(
                    {"row_num": row_num, "start": s_start, "end": s_end}
                )
    for hits in hgt_index.values():
        hits.sort(key=lambda hit: (hit["start"], hit["end"]))
    return hgt_index


def load_roary_gene_instances(roary_file: Path, selected_genes):
    selected_genes = set(selected_genes)
    roary_rows = {}
    with roary_file.open(newline="") as handle:
        reader = csv.DictReader(handle)
        genome_cols = [col.strip('"') for col in reader.fieldnames[14:]]
        for row in reader:
            gene = row["Gene"].strip('"')
            if gene not in selected_genes:
                continue
            genomes = {}
            for col in genome_cols:
                raw_value = row.get(col, "")
                if raw_value:
                    genomes[col] = raw_value.strip('"')
            roary_rows[gene] = genomes
    return roary_rows


def parse_gene_coords(gff_dir: Path, roary_rows, hgt_genomes):
    needed = defaultdict(set)
    for genomes in roary_rows.values():
        for genome, gene_ids in genomes.items():
            if norm_genome(genome) not in hgt_genomes:
                continue
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
        with gff.open() as handle:
            for line in handle:
                if line.startswith("##FASTA"):
                    break
                if line.startswith("#") or "\tCDS\t" not in line or "ID=" not in line:
                    continue
                parts = line.rstrip("\n").split("\t", 8)
                gene_id = gff_attr_id(parts[8])
                if gene_id not in ids:
                    continue
                start, end = sort_interval(parts[3], parts[4])
                coords[(genome, gene_id)] = {
                    "contig": norm_contig(parts[0]),
                    "start": start,
                    "end": end,
                }
    return coords


def hgt_stats_for_gene(gene: str, roary_rows, coords, hgt_index, hgt_genomes):
    gene_instances = set()
    hgt_rows = set()
    max_overlap_bp = 0
    for genome, gene_ids in roary_rows.get(gene, {}).items():
        ngenome = norm_genome(genome)
        if ngenome not in hgt_genomes:
            continue
        for gene_id in gene_ids.split(";"):
            gene_id = gene_id.strip()
            if not gene_id:
                continue
            coord = coords.get((genome, gene_id))
            if not coord:
                continue
            for hit in hgt_index.get((ngenome, coord["contig"]), []):
                if hit["start"] > coord["end"]:
                    break
                if hit["end"] < coord["start"]:
                    continue
                ov = overlap_len(coord["start"], coord["end"], hit["start"], hit["end"])
                if ov <= 0:
                    continue
                gene_instances.add((genome, gene_id))
                hgt_rows.add(hit["row_num"])
                max_overlap_bp = max(max_overlap_bp, ov)
    return len(gene_instances), len(hgt_rows), max_overlap_bp


def direction(beta: float) -> str:
    if beta > 0:
        return "IBD"
    if beta < 0:
        return "Health"
    return "zero_beta"


def analyze_species(species: str, kept_samples: int):
    sig_genes = read_sig_genes(species)
    print(f"  FDR<0.05 genes={len(sig_genes)}", flush=True)
    pangenome_dir = PANGENOME_BASE / species
    roary_file = pangenome_dir / "gene_presence_absence.csv"
    gff_dir = pangenome_dir / "fixed_input_files"
    genome_cols = read_genome_columns(roary_file)
    target_mags = {norm_genome(col) for col in genome_cols}
    hgt_index = load_hgt_index(target_mags)
    print(f"  HGT contigs with hits={len(hgt_index)}", flush=True)
    roary_rows = load_roary_gene_instances(roary_file, {row["variant"] for row in sig_genes})
    hgt_genomes = {genome for genome, _contig in hgt_index}
    coords = parse_gene_coords(gff_dir, roary_rows, hgt_genomes)
    print(f"  gene coordinates={len(coords)}", flush=True)

    detail_rows = []
    summary = {
        "species": species,
        "skani9999_kept_samples": kept_samples,
        "n_fdr05_genes": 0,
        "n_health_sig": 0,
        "n_ibd_sig": 0,
        "n_fdr05_in_hgt": 0,
        "n_health_sig_in_hgt": 0,
        "n_ibd_sig_in_hgt": 0,
    }

    for row in sig_genes:
        beta = row["_beta"]
        direct = direction(beta)
        hgt_gene_instances, hgt_rows, max_overlap_bp = hgt_stats_for_gene(
            row["variant"], roary_rows, coords, hgt_index, hgt_genomes
        )
        in_hgt = hgt_gene_instances > 0
        detail_rows.append(
            {
                "species": species,
                "skani9999_kept_samples": kept_samples,
                "gene": row["variant"],
                "af": row["af"],
                "lrt_pvalue": row["lrt-pvalue"],
                "fdr": row["fdr"],
                "beta": row["beta"],
                "direction": direct,
                "in_hgt_region": "TRUE" if in_hgt else "FALSE",
                "hgt_gene_instances": hgt_gene_instances,
                "hgt_rows": hgt_rows,
                "max_overlap_bp": max_overlap_bp,
            }
        )
        summary["n_fdr05_genes"] += 1
        if direct == "Health":
            summary["n_health_sig"] += 1
            if in_hgt:
                summary["n_health_sig_in_hgt"] += 1
        elif direct == "IBD":
            summary["n_ibd_sig"] += 1
            if in_hgt:
                summary["n_ibd_sig_in_hgt"] += 1
        if in_hgt:
            summary["n_fdr05_in_hgt"] += 1

    return detail_rows, summary


def write_tsv(path: Path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def sanity_check_top50(detail_rows):
    by_species = defaultdict(set)
    for row in detail_rows:
        if row["in_hgt_region"] == "TRUE":
            by_species[row["species"]].add(row["gene"])

    checks = []
    for species, _kept in selected_species():
        summary_path = species_result_dir(species) / "top50_FDR_lt_0.05_within_HGT_overlap_summary.tsv"
        if not summary_path.exists() or summary_path.stat().st_size == 0:
            existing = set()
        else:
            with summary_path.open(newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                existing = {row["group"] for row in reader}
        top50_path = species_result_dir(species) / "IBD_COGs_FDR_lt_0.05_top50.tsv"
        with top50_path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            top50 = {row["variant"] for row in reader}
        recalculated = by_species[species] & top50
        checks.append(
            {
                "species": species,
                "existing_top50_hgt_genes": len(existing),
                "recalculated_top50_hgt_genes": len(recalculated),
                "missing_from_recalculated": ",".join(sorted(existing - recalculated)),
                "extra_in_recalculated": ",".join(sorted(recalculated - existing)),
            }
        )
    return checks


def main():
    detail_rows = []
    summary_rows = []
    for species, kept_samples in selected_species():
        print(f"Analyzing {species} kept={kept_samples}", flush=True)
        detail, summary = analyze_species(species, kept_samples)
        detail_rows.extend(detail)
        summary_rows.append(summary)

    detail_fields = [
        "species",
        "skani9999_kept_samples",
        "gene",
        "af",
        "lrt_pvalue",
        "fdr",
        "beta",
        "direction",
        "in_hgt_region",
        "hgt_gene_instances",
        "hgt_rows",
        "max_overlap_bp",
    ]
    summary_fields = [
        "species",
        "skani9999_kept_samples",
        "n_fdr05_genes",
        "n_health_sig",
        "n_ibd_sig",
        "n_fdr05_in_hgt",
        "n_health_sig_in_hgt",
        "n_ibd_sig_in_hgt",
    ]
    check_fields = [
        "species",
        "existing_top50_hgt_genes",
        "recalculated_top50_hgt_genes",
        "missing_from_recalculated",
        "extra_in_recalculated",
    ]
    write_tsv(PYSEER_BASE / "species_fdr05_gene_hgt_detail.tsv", detail_rows, detail_fields)
    write_tsv(PYSEER_BASE / "species_fdr05_gene_hgt_summary.tsv", summary_rows, summary_fields)
    write_tsv(PYSEER_BASE / "species_fdr05_gene_hgt_top50_sanity.tsv", sanity_check_top50(detail_rows), check_fields)


if __name__ == "__main__":
    main()
