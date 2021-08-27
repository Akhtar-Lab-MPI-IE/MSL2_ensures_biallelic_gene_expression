rule filterFragments_allele:
    input:
        "allelic_bams/{sample}.{suffix}.sorted.bam"
    output:
        shortBAM = temp(os.path.join(short_bams, "allele_specific/{sample}.{suffix}.short.bam")),
        metrics = os.path.join(short_bams, "allele_specific/{sample}.{suffix}.short.metrics")
    log: os.path.join(short_bams, "allele_specific/logs/{sample}.{suffix}.filterFragments.log")
    params:
        maxFragmentSize=maxFragmentSize,
        minFragmentSize=minFragmentSize
    threads: 6
    conda: CONDA_SHARED_ENV
    shell: """
        alignmentSieve --bam {input} \
        --outFile {output.shortBAM} -p {threads} \
        --filterMetrics {output.metrics} \
        --maxFragmentLength {params.maxFragmentSize} \
        --minFragmentLength {params.minFragmentSize} \
        2> {log}
        """


rule filterMetricsToHtml_allele:
    input:
        expand(os.path.join(short_bams, "allele_specific/{sample}.{{suffix}}.short.metrics"), sample=samples)
    output:
        QCrep='Filtering_metrics/Filtering_report.{suffix}.html'
    log:
        err="Filtering_metrics/logs/produce_report.{suffix}.err",
        out="Filtering_metrics/logs/produce_report.{suffix}.out"
    conda: CONDA_RMD_ENV
    threads: 1
    script: "../rscripts/ATACseq_QC_report_template.Rmd"


# MACS2 BAMPE fails if there is just one fragment mapped
rule filterCoveragePerScaffolds_allele:
    input:
        bam = os.path.join(short_bams, "allele_specific/{sample}.{suffix}.short.bam")
    output:
        whitelist = os.path.join(short_bams, "allele_specific/{sample}.{suffix}.chrom.whitelist"),
        shortbai = temp(os.path.join(short_bams, "allele_specific/{sample}.{suffix}.short.bam.bai")),
        bam = os.path.join(short_bams, "allele_specific/{sample}.{suffix}.short.cleaned.bam"),
        bai = os.path.join(short_bams, "allele_specific/{sample}.{suffix}.short.cleaned.bam.bai")
    log: os.path.join(short_bams, "allele_specific/logs/{sample}.{suffix}.filterCoveragePerScaffolds.log")
    params:
        count_cutoff = int(fragmentCountThreshold) * 2 # must contain more than 2 reads, i.e. 1 fragment
    threads: 6
    conda: CONDA_SHARED_ENV
    shell: """
        samtools index -@ {threads} {input.bam} 2> {log}
        samtools idxstats {input.bam} | awk -v cutoff={params.count_cutoff} \'$3 > cutoff\' | cut -f 1 > {output.whitelist} 2>> {log}
        samtools view -@ {threads} -bo {output.bam} {input.bam} $(cat {output.whitelist} | paste -sd\' \') 2>> {log}
        samtools index -@ {threads} {output.bam} 2>> {log}

        """

# MACS2 BAMPE filter: samtools view -b -f 2 -F 4 -F 8 -F 256 -F 512 -F 2048
rule callOpenChromatin_allele:
    input:
        os.path.join(short_bams, "allele_specific/{sample}.{suffix}.short.cleaned.bam")
    output:
        peaks = os.path.join(outdir_MACS2, 'allele_specific/{sample}.{suffix}.filtered.short.BAM_peaks.narrowPeak'),
        xls = os.path.join(outdir_MACS2, 'allele_specific/{sample}.{suffix}.filtered.short.BAM_peaks.xls')
    params:
        directory = os.path.join(outdir_MACS2,'allele_specific'),
        genome_size = int(genome_size),
        name='{sample}',
        qval_cutoff=qval,
        nomodel='--nomodel',
        write_bdg='--bdg',
        fileformat='--format BAMPE'
    threads: 6
    log:
        out = os.path.join(outdir_MACS2,'allele_specific', "logs", "callOpenChromatin", "{sample}.{suffix}_macs2.out"),
        err = os.path.join(outdir_MACS2,'allele_specific', "logs", "callOpenChromatin", "{sample}.{suffix}_macs2.err")
    conda: CONDA_ATAC_ENV
    shell: """
        macs2 callpeak --treatment {input} \
            -g {params.genome_size} \
            --name {params.name}.filtered.short.BAM \
            --outdir {params.directory} \
            {params.fileformat} \
            --qvalue {params.qval_cutoff} \
            {params.nomodel} \
            {params.write_bdg} > {log.out} 2> {log.err}
        """


rule tempChromSizes_allele:
    input: genome_index
    output: temp("HMMRATAC/chrom_sizes")
    log: "HMMRATAC/logs/tempChromSizes.log"
    shell: """
        cut -f 1,2 {input} > {output} 2> {log}
        """


# TODO: -q MINMAPQ -Xmx value is currently hard-coded
# Actually uses 2-4 cores, even though there's no option for it!
# Requires PE data
rule HMMRATAC_peaks_allele:
    input:
        "allelic_bams/{sample}.{suffix}.sorted.bam",
        "allelic_bams/{sample}.{suffix}.sorted.bam.bai",
        "HMMRATAC/chrom_sizes"
    output:
        "HMMRATAC/{sample}.{suffix}.log",
        "HMMRATAC/{sample}.{suffix}.model",
        "HMMRATAC/{sample}.{suffix}_peaks.gappedPeak",
        "HMMRATAC/{sample}.{suffix}_summits.bed",
        "HMMRATAC/{sample}.{suffix}_training.bed"
    log: "HMMRATAC/logs/{sample}.{suffix}.HMMRATAC_peaks.log"
    params:
        blacklist = "-e {}".format(blacklist_bed) if blacklist_bed else ""
    conda: CONDA_ATAC_ENV
    threads: 4
    shell: """
        HMMRATAC -Xmx10G -b {input[0]} -i {input[1]} -g {input[2]} {params.blacklist} -o HMMRATAC/{wildcards.sample} 2> {log}
        """

#Genrich requires namesorted bams
rule namesort_bams_allele:
    input:
        bam = short_bams + "allele_specific/{sample}.{suffix}.short.cleaned.bam"
    output:
        bam = temp(short_bams + "allele_specific/{sample}.{suffix}.short.namesorted.bam")
    log:
        short_bams + "allele_specific/logs/{sample}.{suffix}.namesort.err"
    params:
        tempDir = tempDir
    threads: 4
    conda: CONDA_SAMBAMBA_ENV
    shell: """
        TMPDIR={params.tempDir}
        MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX)
        sambamba sort -t {threads} -o {output.bam} --tmpdir=$MYTEMP -n {input.bam} 2> {log}
        rm -rf $MYTEMP
         """


# Requires PE data
# Should be run once per-group!
rule Genrich_peaks_allele:
    input:
        bams=lambda wildcards: expand(short_bams + "allele_specific/{sample}.{suffix}.short.namesorted.bam", sample=genrichDict[wildcards.group])
    output:
        "Genrich/{group}.{suffix}.narrowPeak"
    log: "Genrich/logs/{group}.{suffix}.Genrich_peaks.log"
    params:
        bams = lambda wildcards: ",".join(expand(short_bams + "allele_specific/{sample}.{suffix}.short.namesorted.bam", sample=genrichDict[wildcards.group])),
        blacklist = "-E {}".format(blacklist_bed) if blacklist_bed else ""
    conda: CONDA_ATAC_ENV
    shell: """
        Genrich  -t {params.bams} -o {output} -r {params.blacklist} -j -y 2> {log}
        """

### deepTools bamCompare subtract #######################################################

rule bamCompare_subtract_allele:
    input:
        bam = os.path.join(short_bams, "allele_specific/{sample}.{suffix}.short.cleaned.bam"),
        bai = os.path.join(short_bams, "allele_specific/{sample}.{suffix}.short.cleaned.bam.bai")
    output:
        "deepTools_ATAC/bamCompare/allele_specific/{sample}.{suffix}.filtered.bw"
    params:
        bwBinSize = bwBinSize,
        genome_size = genome_size,
        ignoreForNorm = "--ignoreForNormalization {}".format(ignoreForNormalization) if ignoreForNormalization else "",
        read_extension = "--extendReads" if pairedEnd else "--extendReads {}".format(fragmentLength),
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed else ""
    log:
        out = "deepTools_ATAC/allele_specific/logs/bamCompare.{sample}.{suffix}.filtered.out",
        err = "deepTools_ATAC/allele_specific/logs/bamCompare.{sample}.{suffix}.filtered.out"
    benchmark:
        "deepTools_ATAC/.benchmark/deepTools_ATAC/logs/bamCompare.{sample}.{suffix}.filtered.benchmark"
    threads: 16
    conda: CONDA_SHARED_ENV
    shell: bamcov_cmd

rule report_allele:
    input:
        #expand("HMMRATAC/{sample}.{suffix}_peaks.gappedPeak", sample=samples, suffix=["genome1", "genome2"]),
        expand(os.path.join(outdir_MACS2, "allele_specific/{sample}.{suffix}.filtered.short.BAM_peaks.xls"), sample = samples, suffix=["genome1", "genome2"]),
        expand(os.path.join(outdir_ATACqc,"allele_specific/{sample}.{suffix}.filtered.BAM_peaks.qc.txt"), sample = samples, suffix=["genome1", "genome2"]),
        expand("deepTools_ATAC/bamCompare/allele_specific/{sample}.{suffix}.filtered.bw",sample=samples, suffix=["genome1", "genome2"]),
        expand(os.path.join(short_bams, "allele_specific/{sample}.{suffix}.short.cleaned.bam"),sample=samples, suffix=["genome1", "genome2"]),
        expand(os.path.join(short_bams, "allele_specific/{sample}.{suffix}.short.cleaned.bam.bai"),sample=samples, suffix=["genome1", "genome2"]),
        expand("Filtering_metrics/allele_specific/Filtering_report.{suffix}.html", suffix=["genome1", "genome2"])
    output:
        "allelic_bams/report.allelic.txt"
    run:
        file = open(output[0],"w")
        file.write("everything about allele specific is done")
        file.close()
        
    