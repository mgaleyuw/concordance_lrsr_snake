import pandas as pd

samples = pd.read_csv(config["targetfile"], header=None, index_col=0)

OUTPUTDIR=config["outputdirectory"]

THREADS=config["threads"]

def get_targets(wildcards):
    f = open(config["targetfile"], "r")
    samplepairs = f.read().split("\n")
    f.close()
    final_targets = [] 
    runType="fullWGS"
    compType="simple"
    if config["bedRestrict"]:
        runType="bed"
    if config["overlapVCF"]:
        compType="full"
    for pair in samplepairs:
        targets=pair.split(" ")
        if(len(targets) != 2):
            raise ValueError("Target file must contain 2 columns of file paths plus a pair identifier")
        trunk=f"{OUTPUTDIR}/{targets[0]}_concordance_{runType}_{compType}"
        final_targets+=[f"{trunk}.genotype_concordance_{x}" for x in ["contingency_metrics", "detail_metrics", "summary_metrics"]]
    return final_targets

def get_inputfiles(wildcards):
    vcf1=samples.loc[wildcards.COMPARISON, 1]
    vcf2=samples.loc[wildcards.COMPARISON, 2]
    return {"vcf1":vcf1, "vcf2":vcf2}
