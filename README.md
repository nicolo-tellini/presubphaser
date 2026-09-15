# presubphaser

Homeologue pairing via shared single-copy BUSCO genes.

`presubphaser` takes a BUSCO `full_table.tsv` from a polyploid assembly and counts, for every pair of chromosomes, how many BUSCO genes they have in common. Homeologous chromosomes (the copies inherited from the two (or more) subgenomes of an allopolyploid) share more BUSCOs with each other than with any other chromosome; the resulting matrix makes the pairing visible.

The output is the input that [SubPhaser](https://github.com/zhangrengang/SubPhaser) expects: SubPhaser polarises homeologous groups into subgenomes, it does not infer the groups themselves.

## Method

Each BUSCO annotated as `Complete` or `Duplicated` is reduced to the set of chromosomes it was found on, one entry per BUSCO per chromosome. From this a presence/absence matrix is built, and the number of BUSCOs shared by each pair of chromosomes is obtained as its cross-product.

The matrix is written as a TSV and drawn as a lower-triangular heatmap with the counts printed in the cells. Homeologous pairs stand out as isolated high-count cells off the diagonal.

Assigning the pairs is left to the user: reading them off the heatmap is fast, unambiguous where the signal is clean. Chromosomes involved in translocations or in a whole-genome duplication older than the allopolyploid event produce several comparable counts instead of one, and those need to be investigated separately.

## Requirements

R ≥ 4.0 with:

- `optparse`
- `gtools`
- `pheatmap`

```r
install.packages(c("optparse", "gtools", "pheatmap"))
```

## Input

The `full_table.tsv` produced by BUSCO in genome mode, found in the run directory as `run_<lineage>/full_table.tsv`. Only the first three columns are used: BUSCO id, status, and sequence name.

Run BUSCO on the assembly you want to phase, with a lineage appropriate for the clade:

```bash
busco -i assembly.fa -l eudicots_odb10 -m genome -c 16 -o busco_out
```

## Usage

```bash
Rscript presubphaser.r --input full_table.tsv --output PREFIX [--chrom_regex REGEX]
```

## Arguments

| Flag | Short | Required | Default | Description |
|---|---|---|---|---|
| `--input` | `-i` | yes | — | Path to the BUSCO `full_table.tsv` |
| `--output` | `-o` | yes | — | Output prefix, without extension |
| `--chrom_regex` | `-r` | no | `SUPER` | Regex selecting which sequence names to keep |

`--chrom_regex` is what keeps unplaced contigs out of the matrix. Set it to whatever prefix your
chromosome-level sequences carry — `SUPER` for a Pretext/YaHS-curated assembly, `chr` if you have
renamed them, `.` to keep everything.

## Examples

Chromosomes named `SUPER_1` … `SUPER_18`:

```bash
Rscript presubphaser.r -i run_eudicots_odb10/full_table.tsv -o hom_pairing
```

Chromosomes renamed to `chr1` … `chr18`:

```bash
Rscript presubphaser.r -i run_eudicots_odb10/full_table.tsv -o hom_pairing -r "^chr"
```

## Output

| File | Description |
|---|---|
| `PREFIX.shared_matrix.tsv` | Square matrix of BUSCOs shared between every pair of chromosomes |
| `PREFIX.heatmap.pdf` | Lower-triangular heatmap of the same matrix, counts printed in the cells |

## Notes

- A BUSCO found more than once on the same chromosome is counted once: the matrix measures shared gene
  content between chromosomes, not copy number within them.
- The diagonal holds the number of BUSCOs on each chromosome and is blanked in the heatmap so that the
  colour scale is set by the off-diagonal counts.
- A clean allotetraploid gives one dominant partner per chromosome and a one-to-one pairing. Quartets of
  chromosomes with comparable counts, rather than pairs, point to a translocation or to an older duplication
  layered under the allopolyploid event — worth resolving with synteny before passing anything to SubPhaser.
- Very low counts across the whole matrix usually mean the regex is keeping fragmented scaffolds, or that
  the BUSCO lineage is too distant for the clade.# presubphaser

