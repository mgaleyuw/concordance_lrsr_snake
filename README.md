# Concordance Workflow

Run concordance between two long read SNV VCFs or a long read and short read VCFs using Picard.
Filter to regions of high reproducibility and low noise; sort output VCF and save the categorized variants separated by membership to only VCF1, only VCF2, differering zygosity, the union, and calls that differ because of differences in how multiallelic records are represented.

## Usage

Supply a three-column targets list (`targets.txt`) with pairs of VCF paths separated by a space to be compared and an identifier for the pair, e.g:
```
M1234_sr_lr M1234.sr.phased.vcf.gz M1234.lr.phased.vcf.gz
M2334_sr_lr M2344.sr.phased.vcf.gz M2344.lr.phased.vcf.gz
M1234_identity M1234.sr.phased.vcf.gz M2344.lr.phased.vcf.gz
```

Modify `config/config.yaml` to point `targetfile` at the `targets.txt` file if you have changed the name.

Specify in the `config/config.yaml` if you would like to keep and categorize VCFs of overlapping/non-overlapping variants by setting `overlapVCF` to `true`. This is usually the goal of comparing SR/LR VCFs. If simply confirming identity using short reads or relatedness between two LR VCFs, set `overlapVCF` to `false`.

Specify in `config/config.yaml` whether analysis should be restricted to a target area, e.g. an exome bed file for comparing SR Exome data to LR data, by setting `bedRestrict` to `true`. Set as `false` when comparing two LR VCFs.

Conda environment is included in `workflow/envs`; the first run will be slow as this environment is built.

Run using `snakemake --use-conda --cores {cores}`

Summarize output for many files using the helper script `post_run/summarizeConcordance.sh`

## Installation

Clone this git repository and run in an environment that includes Snakemake >= 8.16.0.

## Examples 

This example uses the small vcfs in `example` and the current values in `config/config.yaml` to compare data from the same individual sequenced using LR ONT and sequenced using SR ONT.

`overlapVCF` has been set to `true` and `bedRestrict` is set to `true`. 

`restriction_analysis_bed` has been set to the included twist comprehensive exome targets for hg38.

Run:

```
snakemake --use-conda --cores 10
```
