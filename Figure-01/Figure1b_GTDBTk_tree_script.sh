#!/bin/bash
#SBATCH --job-name=gtdbtk_all_v3_rooted
#SBATCH --output=/scratch/hb-fu/haoran/HGT_tree/gtdbtk_all_v3_rooted/gtdbtk_all_v3_rooted.out
#SBATCH --error=/scratch/hb-fu/haoran/HGT_tree/gtdbtk_all_v3_rooted/gtdbtk_all_v3_rooted.err
#SBATCH --partition=regularlong
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=60
#SBATCH --mem=300G

set -euo pipefail

source activate gtdbtk-2.1.1

script_dir=/scratch/p312334/project/10--HGT_isolates/1--General_infomation/phylotree/gtdbtk_all_v3_rooted
wd=/scratch/hb-fu/haoran/HGT_tree/gtdbtk_all_v3_rooted

mkdir -p "${wd}"
export GTDBTK_OUTDIR="${wd}"

python "${script_dir}/generate_gtdbtk_inputs.py"

gtdbtk de_novo_wf \
    --genome_dir "${wd}/genomes" \
    --out_dir "${wd}/output_hpc" \
    --skip_gtdb_refs \
    --cpus "${SLURM_CPUS_PER_TASK}" \
    --outgroup_taxon p__Minisyncoccota \
    --bacteria \
    --custom_taxonomy_file "${wd}/meta_gtdbtk.txt" \
    -x fa
