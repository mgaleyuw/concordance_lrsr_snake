rule run_concordance_simple:
    input:
        unpack(get_inputfiles)
    output:
        contingency="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_simple/{COMPARISON}_concordance_{RUNTYPE}_simple.genotype_concordance.contingency_metrics"]),
        detail="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_simple/{COMPARISON}_concordance_{RUNTYPE}_simple.genotype_concordance.detail_metrics"]),
        summary="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_simple/{COMPARISON}_concordance_{RUNTYPE}_simple.genotype_concordance.summary_metrics"]),
        workdir=temp(directory("".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_simple_tmp"])))
    threads: THREADS
    conda: "../envs/concordance.yaml"
    params:
        scriptname="workflow/scripts/concordance_run.sh",
        outputname="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_simple/{COMPARISON}_concordance_{RUNTYPE}_simple"]),
        reference=config["reference"],
        mask=config["maskfile"],
        bed=config["restriction_analysis_bed"],
        filterPassing=config["filterPassing"]
    log:
        o="logs/{COMPARISON}_{RUNTYPE}_simple_concordance.o.log"
    shell:
        """
        OUTNAME={params.outputname}
        mkdir -p ${{OUTNAME%/*}}
        if [ {wildcards.RUNTYPE} == "restricted" ]
        then
            if [ {params.filterPassing} == "true" ]
            then
                bash {params.scriptname} -1 {input.vcf1} -2 {input.vcf2} -o {params.outputname} -w {output.workdir} -m {params.mask}  -b {params.bed} -p > {log.o} 2>&1
            else
                bash {params.scriptname} -1 {input.vcf1} -2 {input.vcf2} -o {params.outputname} -w {output.workdir} -m {params.mask}  -b {params.bed} > {log.o} 2>&1
            fi
        else
            if [ {params.filterPassing} == "true" ]
            then
                bash {params.scriptname} -1 {input.vcf1} -2 {input.vcf2} -o {params.outputname} -w {output.workdir} -m {params.mask} -p > {log.o} 2>&1
            else
                bash {params.scriptname} -1 {input.vcf1} -2 {input.vcf2} -o {params.outputname} -w {output.workdir} -m {params.mask} > {log.o} 2>&1
            fi
        fi
        """

rule run_concordance_detail:
    input:
        unpack(get_inputfiles)
    output:
        contingency="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_detail/{COMPARISON}_concordance_{RUNTYPE}_detail.genotype_concordance.contingency_metrics"]),
        detail="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_detail/{COMPARISON}_concordance_{RUNTYPE}_detail.genotype_concordance.detail_metrics"]),
        summary="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_detail/{COMPARISON}_concordance_{RUNTYPE}_detail.genotype_concordance.summary_metrics"]),
        vcf="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_detail/{COMPARISON}_concordance_{RUNTYPE}_detail.genotype_concordance.vcf.gz"]),
        vcf_index="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_detail/{COMPARISON}_concordance_{RUNTYPE}_detail.genotype_concordance.vcf.gz.tbi"]),
        workdir=temp(directory("".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_detail_tmp"])))
    threads: THREADS
    conda: "../envs/concordance.yaml"
    params:
        scriptname="workflow/scripts/concordance_run.sh",
        outputname="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_detail/{COMPARISON}_concordance_{RUNTYPE}_detail"]),
        reference=config["reference"],
        mask=config["maskfile"],
        bed=config["restriction_analysis_bed"],
        filterPassing=config["filterPassing"]
    log:
        o="logs/{COMPARISON}_{RUNTYPE}_detail_concordance.o.log"
    shell:
        """
        OUTNAME={params.outputname}
        mkdir -p ${{OUTNAME%/*}}
        if [ {wildcards.RUNTYPE} == "restricted" ]
        then
            if [ {params.filterPassing} == "true" ]
            then
                bash {params.scriptname} -1 {input.vcf1} -2 {input.vcf2} -o {params.outputname} -m {params.mask} -w {output.workdir} -b {params.bed} -f -p > {log.o} 2>&1
            else
                bash {params.scriptname} -1 {input.vcf1} -2 {input.vcf2} -o {params.outputname} -m {params.mask} -w {output.workdir} -b {params.bed} -f > {log.o} 2>&1
            fi
        else
            if [ {params.filterPassing} == "true" ]
            then
                bash {params.scriptname} -1 {input.vcf1} -2 {input.vcf2} -o {params.outputname} -m {params.mask} -w {output.workdir} -f -p > {log.o} 2>&1
            else
                bash {params.scriptname} -1 {input.vcf1} -2 {input.vcf2} -o {params.outputname} -m {params.mask} -w {output.workdir} -f > {log.o} 2>&1
            fi
        fi
        """

rule split_contingency:
    input:
        vcf="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_detail/{COMPARISON}_concordance_{RUNTYPE}_detail.genotype_concordance.vcf.gz"]),
    output:
        ucall_vcf="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_detail/{COMPARISON}_concordance_{RUNTYPE}_detail.genotype_concordance.Unique.",COLNAMES[1],".vcf"]),
        uref_vcf="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_detail/{COMPARISON}_concordance_{RUNTYPE}_detail.genotype_concordance.Unique.",COLNAMES[0],".vcf"]),
        overlap_vcf="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_detail/{COMPARISON}_concordance_{RUNTYPE}_detail.genotype_concordance.Overlap.vcf"]),
        multiallelic_vcf="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_detail/{COMPARISON}_concordance_{RUNTYPE}_detail.genotype_concordance.Multiallelic.Mismatch.vcf"]),
        mixedzygo_vcf="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_detail/{COMPARISON}_concordance_{RUNTYPE}_detail.genotype_concordance.Mixed.Zygosity.Mismatch.vcf"])
    threads: 1
    conda: "../envs/concordance.yaml"
    params: 
        script="workflow/scripts/split_contingency.sh",
        prefix="".join([OUTPUTDIR,"/{COMPARISON}_concordance_{RUNTYPE}_detail/{COMPARISON}_concordance_{RUNTYPE}_detail.genotype_concordance"]),
        vcf1=COLNAMES[0],
        vcf2=COLNAMES[1]
    shell:
        """
        bash {params.script} -i {input.vcf} -o {params.prefix} -1 {params.vcf1} -2 {params.vcf2}
        """