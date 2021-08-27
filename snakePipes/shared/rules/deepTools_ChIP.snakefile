### deepTools bamCompare subtract #######################################################

if bigWigType == "subtract" or bigWigType == "both":
    rule bamCompare_subtract:
        input:
            chip_bam = "filtered_bam/{chip_sample}.filtered.bam",
            chip_bai = "filtered_bam/{chip_sample}.filtered.bam.bai",
            control_bam = lambda wildcards: "filtered_bam/"+get_control(wildcards.chip_sample)+".filtered.bam",
            control_bai = lambda wildcards: "filtered_bam/"+get_control(wildcards.chip_sample)+".filtered.bam.bai"
        output:
            "deepTools_ChIP/bamCompare/{chip_sample}.filtered.subtract.{control_name}.bw"
        params:
            bwBinSize = bwBinSize,
            genome_size = genome_size,
            ignoreForNorm = "--ignoreForNormalization {}".format(ignoreForNormalization) if ignoreForNormalization else "",
            read_extension = "--extendReads" if pairedEnd else "--extendReads {}".format(fragmentLength),
            blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed else "",
            scaleFactors = " --scaleFactorsMethod readCount "  #readCount    SES
        log:
            out = "deepTools_ChIP/logs/{chip_sample}.filtered.subtract.{control_name}.out",
            err = "deepTools_ChIP/logs/{chip_sample}.filtered.subtract.{control_name}.err" 
        benchmark:
            "deepTools_ChIP/.benchmark/{chip_sample}.filtered.subtract.{control_name}.benchmark"
        threads: 16
        conda: CONDA_SHARED_ENV
        shell: bamcompare_subtract_cmd

### deepTools bamCompare log2ratio #######################################################
if bigWigType == "log2ratio" or bigWigType == "both":
    rule bamCompare_log2:
        input:
            chip_bam = "filtered_bam/{chip_sample}.filtered.bam",
            chip_bai = "filtered_bam/{chip_sample}.filtered.bam.bai",
            control_bam = lambda wildcards: "filtered_bam/"+get_control(wildcards.chip_sample)+".filtered.bam",
            control_bai = lambda wildcards: "filtered_bam/"+get_control(wildcards.chip_sample)+".filtered.bam.bai"
        output:
            "deepTools_ChIP/bamCompare/{chip_sample}.filtered.log2ratio.over_{control_name}.bw"
        params:
            bwBinSize = bwBinSize,
            ignoreForNorm = "--ignoreForNormalization {}".format(ignoreForNormalization) if ignoreForNormalization else "",
            read_extension = "--extendReads" if pairedEnd else "--extendReads {}".format(fragmentLength),
            blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed else "",
            scaleFactors = " --scaleFactorsMethod readCount "  #readCount    SES
        log:
            out = "deepTools_ChIP/logs/{chip_sample}.filtered.log2ratio.over_{control_name}.out",
            err = "deepTools_ChIP/logs/{chip_sample}.filtered.log2ratio.over_{control_name}.err"
        benchmark:
            "deepTools_ChIP/.benchmark/{chip_sample}.filtered.log2ratio.over_{control_name}.benchmark"
        threads: 16
        conda: CONDA_SHARED_ENV
        shell: bamcompare_log2_cmd


### deepTools plotEnrichment ###################################################

rule plotEnrichment:
    input:
        bams = expand("filtered_bam/{sample}.filtered.bam", sample = all_samples),
        bais = expand("filtered_bam/{sample}.filtered.bam.bai", sample = all_samples)
    output:
        png = "deepTools_ChIP/plotEnrichment/plotEnrichment.gene_features.png",
        tsv = "deepTools_ChIP/plotEnrichment/plotEnrichment.gene_features.tsv",
    params:
        genes_gtf = genes_gtf,
        labels = " ".join(all_samples),
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed else "",
        read_extension = "--extendReads" if pairedEnd else "--extendReads {}".format(fragmentLength)
    log:
        out = "deepTools_ChIP/logs/plotEnrichment.out",
        err = "deepTools_ChIP/logs/plotEnrichment.err"
    benchmark:
        "deepTools_ChIP/.benchmark/plotEnrichment.benchmark"
    threads: 24
    conda: CONDA_SHARED_ENV
    shell: plotEnrich_chip_cmd


### deepTools plotFingerprint (all files) ######################################

rule plotFingerprint:
    input:
        bams = expand("filtered_bam/{sample}.filtered.bam", sample = all_samples),
        bais = expand("filtered_bam/{sample}.filtered.bam.bai", sample = all_samples)
    output:
        metrics = "deepTools_ChIP/plotFingerprint/plotFingerprint.metrics.txt"
    params:
        labels = " --labels " + " ".join(all_samples),
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed else "",
        read_extension = "--extendReads" if pairedEnd else "--extendReads {}".format(fragmentLength),
        png = "--plotFile deepTools_ChIP/plotFingerprint/plotFingerprint.png" if (len(all_samples)<=20)
              else "",
        jsd = "--JSDsample filtered_bam/{}.filtered.bam".format(control_samples[0]) if (len(control_samples)>0)
            else ""
    log:
        out = "deepTools_ChIP/logs/plotFingerprint.out",
        err = "deepTools_ChIP/logs/plotFingerprint.err"
    benchmark:
        "deepTools_ChIP/.benchmark/plotFingerprint.benchmark"
    threads: 24
    conda: CONDA_SHARED_ENV
    shell: plotFingerprint_cmd
