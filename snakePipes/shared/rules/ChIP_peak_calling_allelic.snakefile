#allelic_samples = [w.replace('allele_specific/', '') for w in chip_samples]
#print(chip_samples)

rule link_bam_bai_allele:
    input:
        bam = "allelic_bams/{allelic_sample}.sorted.bam",
        bai = "allelic_bams/{allelic_sample}.sorted.bam.bai"
    output:
        bam_out = "filtered_bam/allele_specific/{allelic_sample}.filtered.bam",
        bai_out = "filtered_bam/allele_specific/{allelic_sample}.filtered.bam.bai"
    shell: """
        ln -s -f ../../{input.bam} {output.bam_out};
        ln -s -f ../../{input.bai} {output.bai_out}
    """

rule link_bigwig_allele:
    input:
        bigwig = "bamCoverage/allele_specific/{allelic_sample}.seq_depth_norm.bw"
    output:
        bigwig = "bamCoverage/allele_specific/{allelic_sample}.filtered.seq_depth_norm.bw",
        bigwig_out = "bamCoverage/{allelic_sample}.filtered.seq_depth_norm.bw"
    params:
        bigwig = "allele_specific/{allelic_sample}.filtered.seq_depth_norm.bw"
    shell: """
        mv {input.bigwig} {output.bigwig} && ln -s {params.bigwig} bamCoverage/
    """

#rule report_allele:
#    input:
#        expand("filtered_bam/allele_specific/{allelic_sample}.filtered.bam", allelic_sample = [y for #y in allelic_samples if '/allele_specific/' not in y]),
#        expand("filtered_bam/allele_specific/{allelic_sample}.filtered.bam.bai", allelic_sample = [y #for y in allelic_samples if '/allele_specific/' not in y]),
#        expand("bamCoverage/allele_specific/{allelic_sample}.filtered.seq_depth_norm.bw", allelic_sam#ple = [y for y in allelic_samples if '/allele_specific/' not in y]),
#    output:
#        "allelic_bams/report.allelic.txt"
#    run:
#        file = open(output[0],"w")
#        file.write("everything about allele specific is done")
#        file.close()
        
