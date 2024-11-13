# Concordance Workflow

Calculate concordance between two long read SNV VCFs or compare long read to short read VCFs using Picard.

Filter to regions of high reproducibility and low noise; sort output VCF and save the categorized variants separated by membership to only VCF1, only VCF2, differering zygosity, the union, and calls that differ because of differences in how multiallelic records are represented.
 
In the case of exome SR data, restrict analysis to only exons specified in capture. One example bed file is included.

If read depth information is included in the `DP` tag in the VCF, ony positions with >= 5 reads are included. 


## Usage

Supply a three-column targets list (`targets.txt`) with pairs of VCF paths separated by a space to be compared and an identifier for the pair, e.g:

```
M1234_sr_lr M1234.sr.phased.vcf.gz M1234.lr.phased.vcf.gz
M2334_sr_lr M2344.sr.phased.vcf.gz M2344.lr.phased.vcf.gz
M1234_identity M1234.sr.phased.vcf.gz M2344.lr.phased.vcf.gz
```

Modify `config/config.yaml` to point `targetfile` at the `targets.txt` file if you have changed the name.


Specify in the `config/config.yaml` if you would like to keep and categorize VCFs of overlapping/non-overlapping variants by setting `overlapVCF` to `true`. This is usually the goal of comparing SR/LR VCFs. If simply confirming identity using short reads or relatedness between two LR VCFs, set `overlapVCF` to `false`.


If a `PASS` filter has been included in the VCF, this can also be used as a filter by setting `filterPassing` to `true` in the config file. 


Analysis can be restricted to autosomes ony if desired by setting `autosomeOnly` to `true`.


Specify in `config/config.yaml` whether analysis should be restricted to a target area, e.g. an exome bed file for comparing SR Exome data to LR data, by setting `bedRestrict` to `true`. Set as `false` when comparing two LR VCFs.


Set better naming outputs by specifying VCF type names in `columnNames` in `config.yaml`. For example if the first column of VCF files are derived from blood and the second column are derived from saliva, you can specify `blood,saliva` to assign these names to downstream VCFs.


Conda environment is included in `workflow/envs`; the first run will be slow as this environment is built.


Run using `snakemake --use-conda --cores {cores}`


Summarize output for many files using the helper script `post_run/summarizeConcordance.sh`


## Installation

Clone this git repository and run in an environment that includes Snakemake >= 8.16.0.

### Requirements

This workflow has been tested in Rocky Linux and Ubuntu environments only. 


Your OS must be able to install and support the following tools through Conda:

   - bedtools=2.31.1

   - picard=3.1.1

   - bcftools=1.19


## Examples 

This example uses the small vcfs in `example` and the current values in `config/config.yaml` to compare data from the same individual sequenced using LR ONT and sequenced using SR ONT.


`overlapVCF` has been set to `true` and `bedRestrict` is set to `false`.


The SR VCF does not include "PASS" as a filter (it's been prefiltered), so the flag `-p` is not included. 


`restriction_analysis_bed` is not used in this example.


Run:

```
snakemake --use-conda --cores 10
```


In the raw output VCF, the `TRUTH` VCF corresponds to VCF1, in the first column of inputs, and the `CALL` VCF corresponds to VCF2, in the second column of inputs. 


The categorized contingency VCFs are named according to the column names defined in `config.yaml`, where the `TRUTH` VCF has been renamed `shortread` and the `CALL` VCF has been renamed `longread`.

