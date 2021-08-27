### deepTools bamCompare log2ratio #######################################################

if bigWigType == "subtract" or bigWigType == "both":
    rule bamCompare_subtract_allelic:
        input:
            chip_bam = "filtered_bam/allele_specific/{chip_sample}.filtered.bam",
            chip_bai = "filtered_bam/allele_specific/{chip_sample}.filtered.bam.bai",
            control_bam = "filtered_bam/allele_specific/{control_name}.filtered.bam",
            control_bai = "filtered_bam/allele_specific/{control_name}.filtered.bam.bai"
        output:
            "deepTools_ChIP/bamCompare/allele_specific/{chip_sample}.filtered.subtract.{control_name}.bw"
        params:
            bwBinSize = bwBinSize,
            genome_size = genome_size,
            ignoreForNorm = "--ignoreForNormalization {}".format(ignoreForNormalization) if ignoreForNormalization else "",
            read_extension = "--extendReads" if pairedEnd else "--extendReads {}".format(fragmentLength),
            blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed else "",
            scaleFactors = " --scaleFactorsMethod readCount "
        log:
            out = "deepTools_ChIP/logs/allele_specific/{chip_sample}.filtered.subtract.{control_name}.out",
            err = "deepTools_ChIP/logs/allele_specific/{chip_sample}.filtered.subtract.{control_name}.err"
        benchmark:
            "deepTools_ChIP/.benchmark/allele_specific/{chip_sample}.filtered.subtract.{control_name}.benchmark"
        threads: 16
        conda: CONDA_SHARED_ENV
        shell: bamcompare_subtract_cmd

### deepTools bamCompare log2ratio #######################################################
if bigWigType == "log2ratio" or bigWigType == "both":
    rule bamCompare_log2_allelic:
        input:
            chip_bam = "filtered_bam/allele_specific/{chip_sample}.filtered.bam",
            chip_bai = "filtered_bam/allele_specific/{chip_sample}.filtered.bam.bai",
            control_bam = "filtered_bam/allele_specific/{control_name}.filtered.bam",
            control_bai = "filtered_bam/allele_specific/{control_name}.filtered.bam.bai",
        output:
            "deepTools_ChIP/bamCompare/allele_specific/{chip_sample}.filtered.log2ratio.over_{control_name}.bw"
        params:
            bwBinSize = bwBinSize,
            ignoreForNorm = "--ignoreForNormalization {}".format(ignoreForNormalization) if ignoreForNormalization else "",
            read_extension = "--extendReads" if pairedEnd else "--extendReads {}".format(fragmentLength),
            blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed else "",
            scaleFactors = " --scaleFactorsMethod readCount "
        log:
            out = "deepTools_ChIP/logs/allele_specific/{chip_sample}.filtered.log2ratio.over_{control_name}.out",
            err = "deepTools_ChIP/logs/allele_specific/{chip_sample}.filtered.log2ratio.over_{control_name}.err"
        benchmark:
            "deepTools_ChIP/.benchmark/allele_specific/{chip_sample}.filtered.log2ratio.over_{control_name}.benchmark"
        threads: 16
        conda: CONDA_SHARED_ENV
        shell: bamcompare_log2_cmd



### deepTools plotEnrichment ###################################################

rule plotEnrichment_allelic:
    input:
        bams = expand("allelic_bams/{sample}.{suffix}.sorted.bam", sample = all_samples, suffix = ['genome1', 'genome2']),
        bais = expand("allelic_bams/{sample}.{suffix}.sorted.bam.bai", sample = all_samples, suffix = ['genome1', 'genome2'])
    output:
        png = "deepTools_ChIP/plotEnrichment/plotEnrichment.gene_features_allelic.png",
        tsv = "deepTools_ChIP/plotEnrichment/plotEnrichment.gene_features_allelic.tsv",
    conda: CONDA_SHARED_ENV
    params:
        genes_gtf = genes_gtf,
        labels = " ".join(expand("{sample}_{suffix}", sample = all_samples, suffix = ['genome1', 'genome2'])),
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                    else "",
        read_extension = "--extendReads" if pairedEnd
                         else "--extendReads " + str(fragmentLength)
    log:
        out = "deepTools_ChIP/logs/plotEnrichment_allelic.out",
        err = "deepTools_ChIP/logs/plotEnrichment_allelic.err"
    benchmark:
        "deepTools_ChIP/.benchmark/plotEnrichment_allelic.benchmark"
    threads: 24
    shell: plotEnrich_chip_cmd


### deepTools plotFingerprint (all files) ######################################

rule plotFingerprint_allelic:
    input:
        bams = expand("allelic_bams/{sample}.{suffix}.sorted.bam", sample = all_samples, suffix = ['genome1', 'genome2']),
        bais = expand("allelic_bams/{sample}.{suffix}.sorted.bam.bai", sample = all_samples, suffix = ['genome1', 'genome2'])
    output:
        metrics = "deepTools_ChIP/plotFingerprint/plotFingerprint.metrics_allelic.txt"
    conda: CONDA_SHARED_ENV
    params:
        labels = " ".join(expand("{sample}_{suffix}", sample = all_samples, suffix = ['genome1', 'genome2'])),
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                    else "",
        read_extension = "--extendReads" if pairedEnd
                         else "--extendReads " + str(fragmentLength),
        png = "--plotFile deepTools_ChIP/plotFingerprint/plotFingerprint_allelic.png" if (len(all_samples)<=20)
              else "",
        jsd = ""
    log:
        out = "deepTools_ChIP/logs/plotFingerprint_allelic.out",
        err = "deepTools_ChIP/logs/plotFingerprint_allelic.err"
    benchmark:
        "deepTools_ChIP/.benchmark/plotFingerprint_allelic.benchmark"
    threads: 24
    shell: plotFingerprint_cmd
