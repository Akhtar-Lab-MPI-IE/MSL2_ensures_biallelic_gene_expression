import os
import sys

if not os.path.exists(config['genome']):
    sys.exit('{} not found, double check config.'.format(config['genome']))
if not os.path.exists(config['gtf']):
    sys.exit('{} not found, double check config.'.format(config['gtf']))
if not os.path.exists(config['SNPdb']):
    sys.exit('{} not found, double check config.'.format(config['SNPdb']))

CHRS = [
    'chr1',
    'chr2',
    'chr3',
    'chr4',
    'chr5',
    'chr6',
    'chr7',
    'chr8',
    'chr9',
    'chr10',
    'chr11',
    'chr12',
    'chr13',
    'chr14',
    'chr15',
    'chr16',
    'chr17',
    'chr18',
    'chr19',
    'chrX',
    'chrY',
    'chrMT'
]

rule all:
    input:
        # 9sCa
        'STAR_9sCa/SAindex',
        'bt2_9sCa/genome.1.bt2',
        # BlCa
        'STAR_BlCa/SAindex',
        'bt2_BlCa/genome.1.bt2'

rule extractfna:
    output:
        temp('oref/genome.fa')
    params:
        ref = config['genome']
    shell:'''
    gunzip -c {params.ref} > {output}
    '''

rule prep_9sCa:
    input:
        'oref/genome.fa'
    output:
        expand(
            '9sCa/129S1_SvImJ_CAST_EiJ_dual_hybrid.based_on_GRCm38_N-masked/{chr}.N-masked.fa',
            chr=CHRS
        )
    params:
        vcf = config['SNPdb']
    shell:'''
    cd 9sCa
    SNPsplit_genome_preparation --vcf_file {params.vcf} --strain 129S1_SvImJ --strain2 CAST_EiJ --dual_hybrid --full_sequence --reference_genome ../oref
    cd ..
    '''


rule star_9sCa:
    input:
        expand(
            '9sCa/129S1_SvImJ_CAST_EiJ_dual_hybrid.based_on_GRCm38_N-masked/{chr}.N-masked.fa',
            chr=CHRS
        )
    output:
        'STAR_9sCa/SAindex'
    threads: 25
    params:
         gtf = config['gtf'],
         odir = 'STAR_9sCa'
    shell:'''
    STAR --runThreadN {threads} \
      --runMode genomeGenerate \
      --genomeDir {params.odir} \
      --genomeFastaFiles {input} \
      --sjdbGTFfile {params.gtf}
    '''
rule bt2_9sCa:
    input:
        expand(
            '9sCa/129S1_SvImJ_CAST_EiJ_dual_hybrid.based_on_GRCm38_N-masked/{chr}.N-masked.fa',
            chr=CHRS
        )
    output:
        'bt2_9sCa/genome.1.bt2'
    threads: 25
    shell:'''
    cat {input} > bt2_9sCa/genome.fa
    bowtie2-build --threads {threads} bt2_9sCa/genome.fa bt2_9sCa/genome
    '''

rule prep_BlCa:
    input:
        'oref/genome.fa'
    output:
        expand(
            'BlCa/CAST_EiJ_N-masked/{chr}.N-masked.fa',
            chr=CHRS
        )
    params:
        vcf = config['SNPdb']
    shell:'''
    cd BlCa
    SNPsplit_genome_preparation --vcf_file {params.vcf} --strain CAST_EiJ --full_sequence --reference_genome ../oref
    cd ..
    '''

rule star_BlCa:
    input:
        expand(
            'BlCa/CAST_EiJ_N-masked/{chr}.N-masked.fa',
            chr=CHRS
        )
    output:
        'STAR_BlCa/SAindex'
    threads: 25
    params:
         gtf = config['gtf'],
         odir = 'STAR_BlCa'
    shell:'''
    STAR --runThreadN {threads} \
      --runMode genomeGenerate \
      --genomeDir {params.odir} \
      --genomeFastaFiles {input} \
      --sjdbGTFfile {params.gtf}
    '''

rule bt2_BlCa:
    input:
        expand(
            'BlCa/CAST_EiJ_N-masked/{chr}.N-masked.fa',
            chr=CHRS
        )
    output:
        'bt2_BlCa/genome.1.bt2'
    threads: 25
    shell:'''
    cat {input} > bt2_BlCa/genome.fa
    bowtie2-build --threads {threads} bt2_BlCa/genome.fa bt2_BlCa/genome
    '''