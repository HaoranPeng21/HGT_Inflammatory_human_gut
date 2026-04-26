#!/usr/bin/env python3
import csv
import math
from pathlib import Path


BASE = Path("/scratch/p312334/project/10--HGT_isolates/3.2--Pangenome")
SUMMARY = BASE / "species_fdr05_gene_hgt_summary.tsv"
OUT = BASE / "species_fdr05_hgt_fisher.tsv"


def log_choose(n, k):
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_pmf(x, total, success_total, draws):
    return math.exp(
        log_choose(success_total, x)
        + log_choose(total - success_total, draws - x)
        - log_choose(total, draws)
    )


def fisher_greater(a, b, c, d):
    """One-sided Fisher exact test for higher HGT proportion in IBD.

    Table:
        IBD     HGT=a, nonHGT=b
        Health  HGT=c, nonHGT=d
    """
    total = a + b + c + d
    hgt_total = a + c
    ibd_total = a + b
    max_x = min(hgt_total, ibd_total)
    return min(1.0, sum(hypergeom_pmf(x, total, hgt_total, ibd_total) for x in range(a, max_x + 1)))


def odds_ratio(a, b, c, d):
    if b == 0 or c == 0:
        if a == 0 or d == 0:
            return "NA"
        return "Inf"
    return f"{(a * d) / (b * c):.6g}"


def bh_fdr(pvals):
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    qvals = [1.0] * n
    prev = 1.0
    for pos, idx in enumerate(reversed(order), start=1):
        rank = n - pos + 1
        q = min(prev, pvals[idx] * n / rank)
        qvals[idx] = q
        prev = q
    return qvals


def main():
    rows = []
    pvals = []
    with SUMMARY.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            ibd_hgt = int(row["n_ibd_sig_in_hgt"])
            health_hgt = int(row["n_health_sig_in_hgt"])
            ibd_total = int(row["n_ibd_sig"])
            health_total = int(row["n_health_sig"])
            ibd_nonhgt = ibd_total - ibd_hgt
            health_nonhgt = health_total - health_hgt
            p = fisher_greater(ibd_hgt, ibd_nonhgt, health_hgt, health_nonhgt)
            pvals.append(p)
            rows.append(
                {
                    "species": row["species"],
                    "skani9999_kept_samples": row["skani9999_kept_samples"],
                    "health_enriched_genes": health_total,
                    "ibd_enriched_genes": ibd_total,
                    "health_hgt_genes": health_hgt,
                    "ibd_hgt_genes": ibd_hgt,
                    "health_nonhgt_genes": health_nonhgt,
                    "ibd_nonhgt_genes": ibd_nonhgt,
                    "health_hgt_fraction": f"{health_hgt / health_total:.6g}" if health_total else "NA",
                    "ibd_hgt_fraction": f"{ibd_hgt / ibd_total:.6g}" if ibd_total else "NA",
                    "odds_ratio_ibd_vs_health_hgt": odds_ratio(ibd_hgt, ibd_nonhgt, health_hgt, health_nonhgt),
                    "fisher_p_greater_ibd_hgt": f"{p:.6g}",
                }
            )

    qvals = bh_fdr(pvals)
    for row, q in zip(rows, qvals):
        row["fisher_bh_fdr"] = f"{q:.6g}"
        try:
            p = float(row["fisher_p_greater_ibd_hgt"])
            row["nominal_p_lt_0.05"] = "TRUE" if p < 0.05 else "FALSE"
        except ValueError:
            row["nominal_p_lt_0.05"] = "FALSE"
        row["bh_fdr_lt_0.05"] = "TRUE" if q < 0.05 else "FALSE"

    fields = [
        "species",
        "skani9999_kept_samples",
        "health_enriched_genes",
        "ibd_enriched_genes",
        "health_hgt_genes",
        "ibd_hgt_genes",
        "health_nonhgt_genes",
        "ibd_nonhgt_genes",
        "health_hgt_fraction",
        "ibd_hgt_fraction",
        "odds_ratio_ibd_vs_health_hgt",
        "fisher_p_greater_ibd_hgt",
        "fisher_bh_fdr",
        "nominal_p_lt_0.05",
        "bh_fdr_lt_0.05",
    ]
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(OUT)


if __name__ == "__main__":
    main()
