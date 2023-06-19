import os

# Set variables.
fqdir = config['fqdir']
Genomes = config['strains']
DualHybridIndex = config['dualHybridIndex']
genomeFa = config['genomeFa']

def infiles(fqdir):
	samples = []
	for file in os.listdir(fqdir):
		if 'R1' in file:
			samples.append(file.replace('_R1.fastq.gz',''))
	return samples

fqfiles = infiles(fqdir)
localrules: shipTheLogs, prepGenome
rule all:
	input:
		# run FastQC
		expand('FASTQC/{sample}_R1_fastqc.zip', sample=fqfiles),
		# run multiQC
		'multiQC/multiqc_report.html',

		# Create bismark index
		'SNPgenome/Bismark_index/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa',
		# Run Bismark mappings and splitting.
		expand('Bismark/{sample}_R1_bismark_bt2_pe.bam', sample=fqfiles),

		# Split bam files.
		expand('SplitBAM/{sample}.split.txt',sample=fqfiles),
		expand('SplitBAM/{sample}_R1_bismark_bt2_pe.{genome}.bam',genome=Genomes, sample=fqfiles),
		# Run Methylation extractions.
		expand('methExtract/{sample}.extract.txt',sample=fqfiles),
		# fullmethExtract
		expand('fullmethExtract/{sample}.done', sample=fqfiles),
		# Sort bam files
		expand('sortBAM/{sample}.{genome}.sort.bam.bai', genome=Genomes, sample=fqfiles),
		# MethylDackel - mbias
		expand('methylDackel/{sample}.{genome}.CHH.mBias.txt', genome=Genomes, sample=fqfiles),
		# MethylDackel - conversionRates
		expand('methylDackel/{sample}.{genome}.convRate', genome=Genomes, sample=fqfiles),
		# GC Bias
		expand('deepTools/{sample}.{genome}.GCbias.txt', genome=Genomes, sample=fqfiles),
		# Genome coverage
		'deepTools/Cov.txt'

rule fastQC:
	input:
		R1 = fqdir + '{sample}_R1.fastq.gz',
		R2 = fqdir + '{sample}_R2.fastq.gz'
	output: 
		'FASTQC/{sample}_R1_fastqc.zip',
		'FASTQC/{sample}_R2_fastqc.zip'
	threads: 5
	shell:'''
	fastqc -t {threads} -o FASTQC {input.R1} {input.R2}
	'''

rule multiqc:
	input: expand('FASTQC/{samples}_R1_fastqc.zip', samples=fqfiles)
	output: 'multiQC/multiqc_report.html'
	threads: 1
	shell:'''
	multiqc -o multiQC FASTQC
	'''

rule bismarkIndex:
	output: 'SNPgenome/Bismark_index/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa'
	threads: 20
	params: 
		bismarkThreads = 8,
		G1 = Genomes[0],
		G2 = Genomes[1],
		maskedFnaDir = config['Nmaskpath']
	shell:'''
	mkdir -p SNPgenome/Bismark_index
	cat {params.maskedFnaDir}/*fa > SNPgenome/Bismark_index/mm10_Nmasked.fa
	faToTwoBit SNPgenome/Bismark_index/mm10_Nmasked.fa SNPgenome/Bismark_index/mm10_Nmasked.2bit
	bismark_genome_preparation --parallel {params.bismarkThreads} SNPgenome/Bismark_index
	'''

rule bismarkMap:
	input: 
		R1 = fqdir + '{sample}_R1.fastq.gz',
		R2 = fqdir + '{sample}_R2.fastq.gz'
	output: 'Bismark/{sample}_R1_bismark_bt2_pe.bam'
	params: 
		bismarkThreads = 8
	threads: 20
	shell:'''
	bismark -p {params.bismarkThreads} --temp_dir {params.tmpdir} SNPgenome/Bismark_index/ -1 {input.R1} -2 {input.R2} -o Bismark 
	'''

rule dedup:
	input:
		'Bismark/{sample}_R1_bismark_bt2_pe.bam'
	output:
		'dedup/{sample}.deduplicated.bam'
	params:
		sample = '{sample}'
	threads: 5
	shell:'''
	deduplicate_bismark -p --output_dir dedup -o {params.sample} {input}
	'''


rule SNPsplit:
	input: 
		'Bismark/{sample}_R1_bismark_bt2_pe.bam'
	params:
		sample = '{sample}',
		G1 = Genomes[0],
		G2 = Genomes[1]
	output: 
		txt = 'SplitBAM/{sample}.split.txt',
		bam1 = 'SplitBAM/{sample}_R1_bismark_bt2_pe.' + Genomes[0] + '.bam',
		bam2 = 'SplitBAM/{sample}_R1_bismark_bt2_pe.' + Genomes[1] + '.bam',
	threads: 10
	shell:'''
	SNPsplit --snp_file SNPgenome/all_SNPs_CAST_EiJ_GRCm38.txt.gz --bisulfite --paired -o SplitBAM/ {input}
	mv SplitBAM/{params.sample}_R1_bismark_bt2_pe.genome1.bam SplitBAM/{params.sample}_R1_bismark_bt2_pe.{params.G1}.bam
	mv SplitBAM/{params.sample}_R1_bismark_bt2_pe.genome2.bam SplitBAM/{params.sample}_R1_bismark_bt2_pe.{params.G2}.bam
	touch {output}
	'''

rule methExtract:
	input:
		'SplitBAM/{sample}.split.txt',
	params:
		sample = '{sample}',
		bismarkThreads = 6, #per '--parallel thread meth extractor uses ~3 threads.'
		genomeDir = os.path.abspath('SNPgenome/Bismark_index'),
		bam1 = 'SplitBAM/{sample}_R1_bismark_bt2_pe.' + Genomes[0] + '.bam',
		bam2 = 'SplitBAM/{sample}_R1_bismark_bt2_pe.' + Genomes[1] + '.bam'
	output:
		'methExtract/{sample}.extract.txt'
	threads: 20
	shell:'''
	bismark_methylation_extractor --genome_folder {params.genomeDir} --cytosine_report --paired --bedGraph --buffer_size 8G --parallel {params.bismarkThreads} -o methExtract/ {params.bam1}
	bismark_methylation_extractor --genome_folder {params.genomeDir} --cytosine_report --paired --bedGraph --buffer_size 8G --parallel {params.bismarkThreads} -o methExtract/ {params.bam2}
	touch {output}
	'''

rule fullmethExtract:
	input:
		'Bismark/{sample}_R1_bismark_bt2_pe.bam'
	params:
		sample = '{sample}',
		bismarkThreads = 6,
		genomeDir = os.path.abspath('SNPgenome/Bismark_index')
	output:
		'fullmethExtract/{sample}.done'
	threads: 20
	shell:'''
	bismark_methylation_extractor --genome_folder {params.genomeDir} --cytosine_report --paired --bedGraph --buffer_size 8G --parallel {params.bismarkThreads} -o fullmethExtract/ {input}
	touch {output}
	'''

rule bamSort:
	input:
		"SplitBAM/{sample}_R1_bismark_bt2_pe.{genome}.bam"
	output:
		"sortBAM/{sample}.{genome}.sort.bam.bai",
		"sortBAM/{sample}.{genome}.sort.bam"
	params:
		"sortBAM/{sample}.{genome}.sort.bam"
	threads: 10
	shell:'''
	samtools sort -@ 10 -m 2000M -T /data/extended/ {input} -o {params}
	samtools index -@ 10 {params}
	'''

rule methylDackel_CHHbias:
	input:
		'sortBAM/{sample}.{genome}.sort.bam'
	output:
		'methylDackel/{sample}.{genome}.CHH.mBias.txt'
	threads: 10
	shell:'''
	mkdir -p methylDackel
	MethylDackel mbias -@ {threads} --CHH --noCpG --noSVG SNPgenome/Bismark_index/mm10_Nmasked.fa {input} > {output}
	'''

rule ConvRates:
	input:
		'methylDackel/{sample}.{genome}.CHH.mBias.txt'
	output:
		'methylDackel/{sample}.{genome}.convRate'
	params:
		'{genome}'
	threads: 1
	shell:'''
	awk '{{if(NR>1) {{M+=$4; UM+=$5}}}}END{{printf("{wildcards.sample}.{params}\\t%f\\n", 100*(1.0-M/(M+UM)))}}' {input} > {output}
	'''

rule GCBias:
	input:
		'sortBAM/{sample}.{genome}.sort.bam'
	output:
		'deepTools/{sample}.{genome}.GCbias.txt'
	params:
		genomeSize = '2652783500', # Set to mouse genome.
		twobitPath = 'SNPgenome/Bismark_index/mm10_Nmasked.2bit'
	threads: 20
	shell:'''
	mkdir -p deepTools
	computeGCBias -b {input} --effectiveGenomeSize {params.genomeSize} -g {params.twobitPath} -l 300 -o {output} -p {threads}
	'''

rule GenomeCov:
	input:
		expand('sortBAM/{sample}.{genome}.sort.bam', genome=Genomes, sample=fqfiles),
	output:
		'deepTools/Cov.txt'
	threads: 20
	shell:'''
	plotCoverage -b {input} -p {threads} -ct 0 -ct 1 -ct 2 -ct 5 -ct 10 -ct 15 -ct 20 -ct 30 -ct 50 --smartLabels --outCoverageMetrics deepTools/Cov.metrics.txt -o deepTools/Cov.png > {output}
	'''

