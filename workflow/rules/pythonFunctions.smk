import pandas as pd

samples = pd.read_table(config["targetfile"], header=None, index_col=0, sep=" ")

OUTPUTDIR=config["outputdirectory"]
COLNAMES=config["columnNames"].split(",")
THREADS=config["threads"]

def get_targets(wildcards):
    f = open(config["targetfile"], "r")
    samplepairs = f.read().split("\n")
    f.close()
    final_targets = []
    extensions=["contingency_metrics", "detail_metrics", "summary_metrics"] 
    runType="WGS"
    compType="simple"
    if config["bedRestrict"]:
        runType="restricted"
    if config["overlapVCF"]:
        compType="detail"
        extensions+=["vcf.gz", "vcf.gz.tbi", f"Unique.{COLNAMES[1]}.vcf", f"Unique.{COLNAMES[0]}.vcf", "Overlap.vcf", "Multiallelic.Mismatch.vcf", "Mixed.Zygosity.Mismatch.vcf"]
    for pair in samplepairs:
        targets=pair.split(" ")
        if(len(targets) != 3):
            raise ValueError("Target file must contain 2 columns of file paths plus a pair identifier")
        trunk=f"{targets[0]}_concordance_{runType}_{compType}"
        final_targets+=[f"{OUTPUTDIR}/{trunk}/{trunk}.genotype_concordance.{x}" for x in extensions]
    return final_targets

def get_inputfiles(wildcards):
    vcf1=samples.loc[wildcards.COMPARISON, 1]
    vcf2=samples.loc[wildcards.COMPARISON, 2]
    return {"vcf1":vcf1, "vcf2":vcf2}
