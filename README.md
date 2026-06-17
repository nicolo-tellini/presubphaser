# presubphaser

Identify homeologous chromosome groups from BUSCO results by computing pairwise shared single-copy gene counts across chromosomes and visualising them as a heatmap.

## Method

For each pair of chromosomes, the script counts how many single-copy BUSCO genes are present on both (Complete or Duplicated status). This produces a symmetric matrix where high values indicate strong homeology. The matrix is displayed as a lower-triangle heatmap. Optionally, chromosomes are grouped via hierarchical clustering (UPGMA) on the shared-BUSCO distance matrix.

## Requirements

R ≥ 4.0 with the following packages:

```r
install.packages(c("optparse", "gtools", "pheatmap"))
```

## Input

The `full_table.tsv` file produced by [BUSCO](https://busco.ezlab.org/) (v5+). It is the per-sequence summary table found inside the BUSCO output directory:

```
busco_output/
└── run_<lineage>/
    └── full_table.tsv   ← this file
```
It is reccomended using latest busco dataset (v12).

## Usage

```bash
Rscript presubphaser.r --input full_table.tsv --output <prefix> [--chrom_regex <regex>] [--clusters <k>]
```

### Arguments

| Flag | Short | Required | Default | Description |
|---|---|---|---|---|
| `--input` | `-i` | yes | — | Path to BUSCO `full_table.tsv` |
| `--output` | `-o` | yes | — | Output prefix (no extension) |
| `--chrom_regex` | `-r` | no | `chr` | Regex to filter sequence names |
| `--clusters` | `-k` | no | — | Number of homeolog groups; if omitted clustering is skipped |

## Examples

Produce only the heatmap and matrix, then inspect the PDF manually to choose the number of groups:

```bash
Rscript presubphaser.r --input full_table.tsv --output output --chrom_regex chr
```

Once you have decided on a number of groups, run again with `--clusters`:

```bash
Rscript presubphaser.r --input full_table.tsv --output eriogonum --chrom_regex chr --clusters 8
```

## Output

| File | Always produced | Description |
|---|---|---|
| `<prefix>.shared_matrix.tsv` | yes | Symmetric matrix of shared BUSCO counts |
| `<prefix>.heatmap.pdf` | yes | Lower-triangle heatmap (original chromosome order) |
| `<prefix>.clusters.tsv` | only with `--clusters` | Table of chromosome → cluster assignments |
| `<prefix>.heatmap_clustered.pdf` | only with `--clusters` | Heatmap reordered and annotated by cluster |
| `<prefix>.summary.txt` | only with `--clusters` | Per-cluster membership and within-group mean shared BUSCOs |

The summary is also printed to stdout.

## real data example 
We use the genome of Conyza bonariensis (GCA_049639985.1). This plant is a hexaploid with 9 clusters of homeologous chromosomes. We assume we are blind regarding ploidy and number of clusters. In this example, we use Hap 1, which contains a collapsed set of chromosomes for each parental genome. How do we group homeologs? For convenience, homeologs are labeled so it is easier to follow the grouping.

Fist, we run BUSCO to get the `full_table.txt`:

```bash
busco -i GCA_049639985.1_ConBo_ref_v01_genomic.fna -m genome -l eudicotyledons_odb12 -o ConBo_busco -c 16
```

A quick look at the `full_table.txt` can give an hint regarding the ploidy of our samples. 
```bash
grep Duplicates full_table.txt | cut -f 1 | sort | uniq -c | cut -f2 | uniq -c 
```

Then, we run `presubphaser.r` with minimum settings. The --chrom_regex corresponds to the tag on the fasta entry
```bash
Rscript presubphaser.r --input full_table.tsv --output output --chrom_regex CM112840
```
## Notes

This script and documentation were developed with the assistance of Claude (Anthropic). All code was reviewed and validated by the author.

- Only sequences whose names match `--chrom_regex` are retained. Adjust the regex to exclude unanchored scaffolds (e.g. `chr`).
- Clustering uses UPGMA (`hclust(..., method = "average")`) on a distance matrix defined as `1 - (shared / max_shared)`.
- Duplicated BUSCO hits are intentionally included: a BUSCO marked Duplicated across two chromosomes is the core signal of homeology. The `unique()` call only collapses redundant hits of the same BUSCO on the same chromosome (multiple alignment positions on a single sequence), so counts are not inflated.
