#print([y for y in samples if '/allele_specific/' not in y])
rule link_bam_bai_allele:
    input:
        bam = "allelic_bams/{sample}.{suffix}.sorted.bam",
        bai = "allelic_bams/{sample}.{suffix}.sorted.bam.bai"
    output:
        bam_out = "filtered_bam/allele_specific/{sample}.{suffix}.filtered.bam",
        bai_out = "filtered_bam/allele_specific/{sample}.{suffix}.filtered.bam.bai",
    shell: """
        ln -s -f ../../{input.bam} {output.bam_out};
        ln -s -f ../../{input.bai} {output.bai_out}
    """

rule link_bigwig_allele:
    input:
        bigwig = "bamCoverage/allele_specific/{sample}.{suffix}.seq_depth_norm.bw"
    output:
        bigwig = "bamCoverage/allele_specific/{sample}.{suffix}.filtered.seq_depth_norm.bw",
        bigwig_out = "bamCoverage/{sample}.{suffix}.filtered.seq_depth_norm.bw"
    params:
        bigwig = "allele_specific/{sample}.{suffix}.filtered.seq_depth_norm.bw"
    shell: """
        mv {input.bigwig} {output.bigwig} && ln -s {params.bigwig} bamCoverage/
    """

rule filterMetricsToHtml_allelic:
    input:
        expand(os.path.join(short_bams, "allele_specific/{sample}.{suffix}.short.metrics"), sample = [y for y in samples if '/allele_specific/' not in y], suffix=["genome1", "genome2"])
    output:
        QCrep='Filtering_metrics/Filtering_report.allelic.html'
    log:
        err="Filtering_metrics/logs/produce_report.allelic.err",
        out="Filtering_metrics/logs/produce_report.allelic.out"
    conda: CONDA_RMD_ENV
    threads: 1
    script: "../rscripts/ATACseq_QC_report_template.Rmd"


def run_HMMRATAC_allelic():
    file_list = []
    if peakCaller == "HMMRATAC":
        file_list = expand("HMMRATAC/allele_specific/{sample}.{suffix}_peaks.gappedPeak", sample = [y for y in samples if '/allele_specific/' not in y], suffix=["genome1", "genome2"])
    return (file_list)

def run_genrich_allelic():
    if peakCaller == "Genrich":
        file_list = ["Genrich/all_samples.narrowPeak"]
        if sampleSheet:
            file_list = ["Genrich/{}.narrowPeak".format(x) for x in genrichDict.keys()]
        file_list.append(expand(short_bams + "allele_specific/{sample}.{suffix}.short.namesorted.bam",sample = [y for y in samples if '/allele_specific/' not in y], suffix=["genome1", "genome2"]))
        return (file_list)
    return ([])

def run_ATAC_allelic():
    file_list = []
    if os.path.isdir('allelic_bams') and os.listdir('allelic_bams') != []:
        file_list.append( expand("filtered_bam/allele_specific/{sample}.{suffix}.filtered.bam", sample = [y for y in samples if '/allele_specific/' not in y], suffix=["genome1", "genome2"]) )
        file_list.append( expand("filtered_bam/allele_specific/{sample}.{suffix}.filtered.bam.bai", sample = [y for y in samples if '/allele_specific/' not in y], suffix=["genome1", "genome2"]) )
        file_list.append( expand("bamCoverage/allele_specific/{sample}.{suffix}.filtered.seq_depth_norm.bw", sample = [y for y in samples if '/allele_specific/' not in y], suffix=["genome1", "genome2"]) )
        file_list.append( expand(os.path.join(outdir_MACS2, "allele_specific/{sample}.{suffix}.filtered.short.BAM_peaks.xls"), sample = [y for y in samples if '/allele_specific/' not in y], suffix=["genome1", "genome2"]) )
        file_list.append( expand(os.path.join(outdir_ATACqc,"allele_specific/{sample}.{suffix}.filtered.BAM_peaks.qc.txt"), sample = [y for y in samples if '/allele_specific/' not in y], suffix=["genome1", "genome2"]) )
        file_list.append( expand("deepTools_ATAC/bamCompare/allele_specific/{sample}.{suffix}.filtered.bw",sample = [y for y in samples if '/allele_specific/' not in y], suffix=["genome1", "genome2"]) )
        file_list.append( expand(os.path.join(short_bams, "allele_specific/{sample}.{suffix}.short.cleaned.bam"),sample = [y for y in samples if '/allele_specific/' not in y], suffix=["genome1", "genome2"]) )
        file_list.append( expand(os.path.join(short_bams, "allele_specific/{sample}.{suffix}.short.cleaned.bam.bai"),sample = [y for y in samples if '/allele_specific/' not in y], suffix=["genome1", "genome2"]) )
        file_list.append( 'Filtering_metrics/Filtering_report.allelic.html' )
    return(file_list)


#rule report_allele:
#    input:
#        ## run csaw if asked for
#        #run_CSAW_allelic(),
#    output:
#        "allelic_bams/report.allelic.txt"
#    run:
#        file = open(output[0],"w")
#        file.write("everything about allele specific is done")
#        file.close()
        
    