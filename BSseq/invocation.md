### setup and run BSseq data.

create conda environment:

  > conda env create -f condaEnv.yml
  > conda activate alleleWGBS

Set up the config (snakemakeConfig.yaml):

  - fqdir: path to fastq files.
  - strains: strains included. If one of the strains is C57BL_6NJ it should be the first in the list.
  - dualHybridIndex: boolean to indicate dual hybrid mode (False in case none of the strains is C57BL_6NJ).
  - genomeFa: path to mouse reference genome (fasta format).

Invocation:

 > snakemake -s Allele.Snakefile --configfile snakemakeConfig.yaml -c 5 --latency-wait 300

Output:


