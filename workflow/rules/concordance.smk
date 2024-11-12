rule run_concordance_simple:
    input:
        unpack(get_inputfiles)
    output:
        contingency="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_simple.genotype_concordance.contingency_metrics"]),
        detail="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_simple.genotype_concordance.detail_metrics"]),
        summary="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_simple.genotype_concordance.summary_metrics"]),
        workdir=temp(directory("".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_simple_tmp"])))
    threads: THREADS
    conda: "../envs/concordance.yaml"
    params:
        scriptname=branch(evaluate("{RUNTYPE}=='fullWGS'"), then="workflow/scripts/concordance_run.sh", otherwise="workflow/scripts/concordance_run_bed.sh"),
        outputname="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_simple"]),
        reference=config["reference"],
        mask=config["maskfile"],
        bed=config["restriction_analysis_bed"]
    log:
        o="logs/{COMPARISON}_{RUNTYPE}_simple_concordance.o.log"
    shell:
        """
        if [ {wildcards.RUNTYPE} == "bed" ]
        then
            bash {params.scriptname} -1 {input.vcf1} -2 {input.vcf2} -o {params.outputname} -m {params.mask}  -b {params.bed} > {log.o} 2>&1
        else
            bash {params.scriptname} -1 {input.vcf1} -2 {input.vcf2} -o {params.outputname} -m {params.mask} > {log.o} 2>&1
        fi
        """

rule run_concordance_full:
    input:
        unpack(get_inputfiles)
    output:
        contingency="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_full.genotype_concordance.contingency_metrics"]),
        detail="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_full.genotype_concordance.detail_metrics"]),
        summary="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_full.genotype_concordance.summary_metrics"]),
        vcf="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_full.genotype_concordance.vcf.gz"]),
        vcf_index="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_full.genotype_concordance.vcf.gz.tbi"]),
        workdir=temp(directory("".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_full_tmp"])))
    threads: THREADS
    conda: "../envs/concordance.yaml"
    params:
        scriptname=branch(evaluate("{RUNTYPE}=='fullWGS'"), then="workflow/scripts/concordance_run.sh", otherwise="workflow/scripts/concordance_run_bed.sh"),
        outputname="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_full"]),
        reference=config["reference"],
        mask=config["maskfile"],
        bed=config["restriction_analysis_bed"]
    log:
        o="logs/{COMPARISON}_{RUNTYPE}_full_concordance.o.log"
    shell:
        """
        if [ {wildcards.RUNTYPE} == "bed" ]
        then
            bash {params.scriptname} -1 {input.vcf1} -2 {input.vcf2} -o {params.outputname} -m {params.mask}  -b {params.bed} -f > {log.o} 2>&1
        else
            bash {params.scriptname} -1 {input.vcf1} -2 {input.vcf2} -o {params.outputname} -m {params.mask} -f > {log.o} 2>&1
        fi
        """

rule split_contingency:
    input:
        vcf="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_full.genotype_concordance.vcf.gz"]),
    output:
        ucall_vcf="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_full.genotype_concordance.Unique.Call.vcf"]),
        uref_vcf="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_full.genotype_concordance.Unique.Reference.vcf"]),
        overlap_vcf="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_full.genotype_concordance.Overlap.Call.Reference.vcf"]),
        multiallelic_vcf="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_full.genotype_concordance.Multiallelic.Mismatch.vcf"]),
        mixedzygo_vcf="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_full.genotype_concordance.Mixed.Zygosity.Mismatch.vcf"])
    threads: 1
    conda: "../envs/concordance.yaml"
    params: 
        script="workflow/split_contingency.sh",
        prefix="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_full.genotype_concordance"])
    shell:
        """
        bash {params.script} -i {input.vcf} -o {params.prefix}
        """