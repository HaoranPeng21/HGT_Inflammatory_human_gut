Figure 1c. Overlap of horizontally transferred ARG, CAZyme, and virulence-associated functions across species pairs.

This panel is an UpSet plot summarizing whether each species pair had horizontally transferred sequences annotated with any of three functional classes: antimicrobial resistance genes (ARG), carbohydrate-active enzymes (CAZyme), and virulence-associated genes. The analysis was performed at the species-pair level. Genome-pair-level HGT annotations were first read from `data/all_within_individual_cross_species_genome_pairs_2000bp_hgt_flag_plasmid_func.csv`, which contains 546,926 genome-pair rows and 8,217 unique species pairs. For each species pair, each functional class was converted to a binary presence/absence indicator by taking the maximum across all genome pairs belonging to that species pair.

Functional-class definitions were based on non-empty annotation fields in the genome-pair table. ARG presence was defined by a non-empty `Resfams_core_v1.2` annotation, CAZyme presence by a non-empty `CAZyme` annotation, and virulence presence by a non-empty `VFDB_setB_20260327` annotation. Empty strings, `nan`, and `None` were treated as absence. A species pair was counted as present for a class if at least one genome pair within that species pair contained an HGT annotation in that class.

The top bar plot shows the number of species pairs in each exact intersection of functional classes. The x-axis columns correspond to the intersections shown in the dot matrix below; black dots indicate functional classes included in that intersection, grey dots indicate absent classes, and black vertical lines connect functional classes that co-occur in the same intersection. The y-axis of the top bar plot is the species-pair count. The horizontal bar plot on the right shows the total set size for each functional class, i.e. the number of species pairs containing at least one HGT annotation in that class, regardless of overlap with the other classes. The plotted UpSet intersections exclude the `None` group, which contained species pairs with no ARG, CAZyme, or virulence annotation.

Across 8,217 species pairs, 262 had at least one of the three functional classes. The total set sizes were 195 ARG-positive species pairs, 151 CAZyme-positive species pairs, and 109 virulence-positive species pairs. Exact intersections were: ARG only, 87 species pairs; ARG + CAZyme + Virulence, 67; CAZyme only, 34; ARG + CAZyme, 32; CAZyme + Virulence, 18; Virulence only, 15; and ARG + Virulence, 9. The remaining 7,955 species pairs had none of the three functional classes and were not drawn as an UpSet column.

No statistical test is shown in this panel. The plot is a descriptive overlap summary of functional annotations among HGT-positive species-pair comparisons.

Source files archived for this panel:

- Plotted intersection-count table: `README/Figure1c_Table.csv`
- Species-pair presence/absence table: `README/Figure1c_Table_presence.csv`
- Summary of counts: `README/Figure1c_Table_summary.tsv`
- Plotting script: `README/Figure1c_script.py`
- Species-pair table-building script: `README/Figure1c_build_table_script.py`
- Genome-pair functional annotation aggregation script: `README/Figure1c_genome_pair_annotation_script.py`
