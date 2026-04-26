import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

fix = Path('/scratch/p312334/project/10--HGT_isolates/1--General_infomation/Fixation_100')
outdir = fix / 'directional_prevalence_v4.1_skani9999'
summary = Path('/scratch/p312334/project/10--HGT_isolates/data/__binstoget_isolates_summary_v4.1_skani9999.csv')
gene_pairs = fix / 'HGT_table_2000bp_isolates_within_plasmid.gene_pairs_with_individual_single.filtered_by_summary_v4.1_skani9999.csv'
pan_root = Path('/scratch/p312334/project/10--HGT_isolates/pangenome_100')

individual_out = outdir / 'HGT_table_2000bp_isolates_within_plasmid.gene_individual_prevalence.directional.csv'
side_out = outdir / 'HGT_table_2000bp_isolates_within_plasmid.gene_prevalence_by_side_cluster_species.csv'
pair_out = outdir / 'HGT_table_2000bp_isolates_within_plasmid.gene_pair_prevalence_by_directional_cluster_species.csv'

outdir.mkdir(parents=True, exist_ok=True)

ind_species_genomes = defaultdict(set)
with summary.open(newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
        ind_species_genomes[(row['MGS_ID'], row['secondary_cluster'])].add(row['Genome_file'] + '.unicycler')

species_dirs = {}
for d in pan_root.iterdir():
    if d.is_dir():
        m = re.search(r'__(\d+_\d+)$', d.name)
        if m:
            species_dirs[m.group(1)] = d

needed = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
pair_counts = Counter()
rows = skipped_side_empty = skipped_pair_empty = 0
with gene_pairs.open(newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
        rows += 1
        individual = row['individual']
        c1 = row['cluster1']
        s1 = row['species_id1']
        c2 = row['cluster2']
        s2 = row['species_id2']

        if c1:
            needed['1'][s1][c1].add(individual)
        else:
            skipped_side_empty += 1
        if c2:
            needed['2'][s2][c2].add(individual)
        else:
            skipped_side_empty += 1

        if c1 and c2:
            pair_counts[(c1, s1, c2, s2)] += 1
        else:
            skipped_pair_empty += 1

print(f'input gene-pair rows: {rows}')
print(f'skipped side entries with empty cluster: {skipped_side_empty}')
print(f'unique directional pair keys: {len(pair_counts)}')
print(f'skipped pair rows with empty cluster1 or cluster2: {skipped_pair_empty}')

needed_by_species = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
for side, species_to_clusters in needed.items():
    for species_id, cluster_to_inds in species_to_clusters.items():
        for cluster, individuals in cluster_to_inds.items():
            needed_by_species[species_id][cluster][side].update(individuals)

fieldnames_ind = [
    'side', 'cluster', 'species_id', 'individual', 'gene_genome_count',
    'total_genome_count', 'individual_gene_prev', 'summary_total_genome_count',
    'missing_genome_count_in_pan'
]

written = missing_species_dir = missing_cluster = 0
with individual_out.open('w', newline='') as fout:
    writer = csv.DictWriter(fout, fieldnames=fieldnames_ind)
    writer.writeheader()
    for si, (species_id, cluster_to_sides) in enumerate(sorted(needed_by_species.items()), 1):
        species_dir = species_dirs.get(species_id)
        n_requests = sum(len(individuals) for side_to_inds in cluster_to_sides.values() for individuals in side_to_inds.values())
        if species_dir is None:
            missing_species_dir += n_requests
            continue
        gpa = species_dir / 'gene_presence_absence.csv'
        if not gpa.exists():
            missing_species_dir += n_requests
            continue

        requested_genomes = set()
        for side_to_inds in cluster_to_sides.values():
            for individuals in side_to_inds.values():
                for individual in individuals:
                    requested_genomes.update(ind_species_genomes.get((individual, species_id), set()))

        with gpa.open(newline='') as f:
            reader = csv.reader(f)
            header = next(reader)
            header_index = {name: idx for idx, name in enumerate(header)}
            genome_to_col = {g: header_index[g] for g in requested_genomes if g in header_index}
            needed_clusters = set(cluster_to_sides)
            cluster_rows = {}
            for row in reader:
                if row and row[0] in needed_clusters:
                    cluster_rows[row[0]] = row
                    if len(cluster_rows) == len(needed_clusters):
                        break

        for cluster, side_to_inds in sorted(cluster_to_sides.items()):
            gpa_row = cluster_rows.get(cluster)
            if gpa_row is None:
                missing_cluster += sum(len(individuals) for individuals in side_to_inds.values())
                continue
            for side, individuals in sorted(side_to_inds.items()):
                for individual in sorted(individuals):
                    summary_genomes = ind_species_genomes.get((individual, species_id), set())
                    pan_genomes = [g for g in summary_genomes if g in genome_to_col]
                    total = len(pan_genomes)
                    present = 0
                    for genome in pan_genomes:
                        idx = genome_to_col[genome]
                        if idx < len(gpa_row) and gpa_row[idx].strip():
                            present += 1
                    prev = present / total if total else ''
                    writer.writerow({
                        'side': side,
                        'cluster': cluster,
                        'species_id': species_id,
                        'individual': individual,
                        'gene_genome_count': present,
                        'total_genome_count': total,
                        'individual_gene_prev': f'{prev:.6g}' if isinstance(prev, float) else '',
                        'summary_total_genome_count': len(summary_genomes),
                        'missing_genome_count_in_pan': len(summary_genomes) - total,
                    })
                    written += 1
        print(f'processed species {si}/{len(needed_by_species)}: {species_id}; written {written}', flush=True)

print(f'wrote individual directional prevalence: {individual_out}')
print(f'individual directional prevalence rows: {written}')
print(f'missing species-dir requests: {missing_species_dir}')
print(f'missing cluster requests: {missing_cluster}')

agg = defaultdict(lambda: [0, 0, 0])
with individual_out.open(newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
        key = (row['side'], row['cluster'], row['species_id'])
        agg[key][0] += int(row['gene_genome_count'])
        agg[key][1] += int(row['total_genome_count'])
        agg[key][2] += 1

with side_out.open('w', newline='') as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            'side', 'cluster', 'species_id', 'gene_genome_count',
            'total_genome_count', 'gene_prevalence', 'individual_count'
        ],
    )
    writer.writeheader()
    for (side, cluster, species_id), (gene_count, total_count, individual_count) in sorted(agg.items()):
        prev = gene_count / total_count if total_count else ''
        writer.writerow({
            'side': side,
            'cluster': cluster,
            'species_id': species_id,
            'gene_genome_count': gene_count,
            'total_genome_count': total_count,
            'gene_prevalence': f'{prev:.6g}' if isinstance(prev, float) else '',
            'individual_count': individual_count,
        })
print(f'wrote side-level prevalence: {side_out}')
print(f'side-level prevalence rows: {len(agg)}')

with pair_out.open('w', newline='') as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            'cluster1', 'species_id1', 'gene1_genome_count', 'gene1_total_genome_count',
            'gene1_prevalence', 'gene1_individual_count',
            'cluster2', 'species_id2', 'gene2_genome_count', 'gene2_total_genome_count',
            'gene2_prevalence', 'gene2_individual_count', 'gene_pair_row_count'
        ],
    )
    writer.writeheader()
    for (c1, s1, c2, s2), count in sorted(pair_counts.items()):
        a1 = agg.get(('1', c1, s1), [0, 0, 0])
        a2 = agg.get(('2', c2, s2), [0, 0, 0])
        p1 = a1[0] / a1[1] if a1[1] else ''
        p2 = a2[0] / a2[1] if a2[1] else ''
        writer.writerow({
            'cluster1': c1,
            'species_id1': s1,
            'gene1_genome_count': a1[0],
            'gene1_total_genome_count': a1[1],
            'gene1_prevalence': f'{p1:.6g}' if isinstance(p1, float) else '',
            'gene1_individual_count': a1[2],
            'cluster2': c2,
            'species_id2': s2,
            'gene2_genome_count': a2[0],
            'gene2_total_genome_count': a2[1],
            'gene2_prevalence': f'{p2:.6g}' if isinstance(p2, float) else '',
            'gene2_individual_count': a2[2],
            'gene_pair_row_count': count,
        })
print(f'wrote directional pair prevalence: {pair_out}')
print(f'directional pair prevalence rows: {len(pair_counts)}')
