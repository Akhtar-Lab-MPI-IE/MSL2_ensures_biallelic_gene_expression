import os

RNASAMPLES, = glob_wildcards(config['fqdir_rna'] + '/{sample}_R1.fastq.gz')
ATACSAMPLES, = glob_wildcards(config['fqdir_atac'] + '/{sample}_a_I1.fastq.gz')
IX=['a', 'b', 'c', 'd']

rule all:
  input:
    expand(
      'splitFQ_RNA/{rnasample}_g1_S1_L001_R1_001.fastq.gz',
      rnasample=RNASAMPLES
    ),
    expand(
      'splitFQ_ATAC/{atacsample}_{ix}_g1_S1_L001_R1_001.fastq.gz',
      atacsample=ATACSAMPLES,
      ix=IX
    )

rule STAR:
  input:
    config['fqdir_rna'] + '/{rnasample}_R2.fastq.gz'
  output:
    temp('STAR/{rnasample}.bam')
  params:
    gdir = config['star_ix'],
    sname = lambda wildcards: 'STAR/' + wildcards.rnasample
  threads: 25
  shell:'''
  STAR --genomeDir {params.gdir} \
    --runThreadN {threads} \
    --readFilesCommand zcat \
    --outFileNamePrefix {params.sname} \
    --readFilesIn {input} \
    --outSAMattributes NH HI NM MD \
    --outSAMtype BAM Unsorted
  mv {params.sname}Aligned.out.bam {params.sname}.bam
  '''

rule bt2:
  input:
    R1 = config['fqdir_atac'] + '/{atacsample}_{ix}_R1.fastq.gz',
    R2 = config['fqdir_atac'] + '/{atacsample}_{ix}_R2.fastq.gz'
  output:
    temp('bowtie2/{atacsample}_{ix}.bam')
  params:
    btix = config['bt_ix']
  threads: 25
  shell:'''
  bowtie2 -X 1000 -x {params.btix} -1 {input.R1} -2 {input.R2} -p {threads} | samtools view -Sb - | samtools sort -n -m 25G -@ 2 -O bam - > {output}
  '''

rule splitATAC:
  input:
    'bowtie2/{atacsample}_{ix}.bam'
  output:
    g1 = temp('SNPsplit_ATAC/{atacsample}_{ix}.genome1.bam'),
    g2 = temp('SNPsplit_ATAC/{atacsample}_{ix}.genome2.bam')
  params:
    snpfile = config['snpfile']
  threads: 15
  shell:'''
  SNPsplit --snp_file {params.snpfile} --no_sort --paired -o SNPsplit_ATAC {input}
  '''

rule splitRNA:
  input:
    'STAR/{rnasample}.bam'
  output:
    g1 = temp('SNPsplit_RNA/{rnasample}.genome1.bam'),
    g2 = temp('SNPsplit_RNA/{rnasample}.genome2.bam')
  params:
    snpfile = config['snpfile']
  threads: 15
  shell:'''
  SNPsplit --single_end --snp_file {params.snpfile} -o SNPsplit_RNA {input}
  '''

rule ixATAC:
  input:
    g1 = 'SNPsplit_ATAC/{atacsample}_{ix}.genome1.bam',
    g2 = 'SNPsplit_ATAC/{atacsample}_{ix}.genome2.bam'
  output:
    g1 = temp('SNPsplit_ATAC/{atacsample}_{ix}.g1.ix.gz'),
    g2 = temp('SNPsplit_ATAC/{atacsample}_{ix}.g2.ix.gz')
  threads: 2
  shell:'''
  samtools view {input.g1} | cut -f1 | gzip > {output.g1}
  samtools view {input.g2} | cut -f1 | gzip > {output.g2}
  '''

rule ixRNA:
  input:
    g1 = 'SNPsplit_RNA/{rnasample}.genome1.bam',
    g2 = 'SNPsplit_RNA/{rnasample}.genome2.bam'
  output:
    g1 = temp('SNPsplit_RNA/{rnasample}.g1.ix.gz'),
    g2 = temp('SNPsplit_RNA/{rnasample}.g2.ix.gz')
  threads: 2
  shell:'''
  samtools view {input.g1} | cut -f1 | gzip > {output.g1}
  samtools view {input.g2} | cut -f1 | gzip > {output.g2}
  '''

rule splitfq_rna:
  input:
    g1 = 'SNPsplit_RNA/{rnasample}.g1.ix.gz',
    g2 = 'SNPsplit_RNA/{rnasample}.g2.ix.gz'
  output:
    R1_g1 = 'splitFQ_RNA/{rnasample}_g1_S1_L001_R1_001.fastq.gz',
    R2_g1 = 'splitFQ_RNA/{rnasample}_g1_S1_L001_R2_001.fastq.gz',
    R1_g2 = 'splitFQ_RNA/{rnasample}_g2_S1_L001_R1_001.fastq.gz',
    R2_g2 = 'splitFQ_RNA/{rnasample}_g2_S1_L001_R2_001.fastq.gz'
  threads: 8
  params:
    spfq = config['splitfq'],
    R1 = lambda wildcards: config['fqdir_rna'] + '/{}_R1.fastq.gz'.format(wildcards.rnasample),
    R2 = lambda wildcards: config['fqdir_rna'] + '/{}_R2.fastq.gz'.format(wildcards.rnasample)
  shell:'''
  {params.spfq} --fq {params.R1} --ix {input.g1} | gzip > {output.R1_g1}&
  {params.spfq} --fq {params.R2} --ix {input.g1} | gzip > {output.R2_g1}&
  {params.spfq} --fq {params.R1} --ix {input.g2} | gzip > {output.R1_g2}&
  {params.spfq} --fq {params.R2} --ix {input.g2} | gzip > {output.R2_g2}&
  '''

rule splitfq_atac:
  input:
    g1 = 'SNPsplit_ATAC/{atacsample}_{ix}.g1.ix.gz',
    g2 = 'SNPsplit_ATAC/{atacsample}_{ix}.g2.ix.gz'
  output:
    R1_g1 = 'splitFQ_ATAC/{atacsample}_{ix}_g1_S1_L001_R1_001.fastq.gz',
    R2_g1 = 'splitFQ_ATAC/{atacsample}_{ix}_g1_S1_L001_R2_001.fastq.gz',
    I1_g1 = 'splitFQ_ATAC/{atacsample}_{ix}_g1_S1_L001_I1_001.fastq.gz',
    I2_g1 = 'splitFQ_ATAC/{atacsample}_{ix}_g1_S1_L001_I2_001.fastq.gz',
    R1_g2 = 'splitFQ_ATAC/{atacsample}_{ix}_g2_S1_L001_R1_001.fastq.gz',
    R2_g2 = 'splitFQ_ATAC/{atacsample}_{ix}_g2_S1_L001_R2_001.fastq.gz',
    I1_g2 = 'splitFQ_ATAC/{atacsample}_{ix}_g2_S1_L001_I1_001.fastq.gz',
    I2_g2 = 'splitFQ_ATAC/{atacsample}_{ix}_g2_S1_L001_I2_001.fastq.gz'
  threads: 16
  params:
    spfq = config['splitfq'],
    R1 = lambda wildcards: config['fqdir_atac'] + '/{}_{}_R1.fastq.gz'.format(wildcards.atacsample, wildcards.ix),
    R2 = lambda wildcards: config['fqdir_atac'] + '/{}_{}_R2.fastq.gz'.format(wildcards.atacsample, wildcards.ix),
    I1 = lambda wildcards: config['fqdir_atac'] + '/{}_{}_I1.fastq.gz'.format(wildcards.atacsample, wildcards.ix),
    I2 = lambda wildcards: config['fqdir_atac'] + '/{}_{}_I2.fastq.gz'.format(wildcards.atacsample, wildcards.ix)
  shell:'''
  {params.spfq} --fq {params.R1} --ix {input.g1} | gzip > {output.R1_g1}&
  {params.spfq} --fq {params.R2} --ix {input.g1} | gzip > {output.R2_g1}&
  {params.spfq} --fq {params.I1} --ix {input.g1} | gzip > {output.I1_g1}&
  {params.spfq} --fq {params.I2} --ix {input.g1} | gzip > {output.I2_g1}&
  {params.spfq} --fq {params.R1} --ix {input.g2} | gzip > {output.R1_g2}&
  {params.spfq} --fq {params.R2} --ix {input.g2} | gzip > {output.R2_g2}&
  {params.spfq} --fq {params.I1} --ix {input.g2} | gzip > {output.I1_g2}&
  {params.spfq} --fq {params.I2} --ix {input.g2} | gzip > {output.I2_g2}
  '''
