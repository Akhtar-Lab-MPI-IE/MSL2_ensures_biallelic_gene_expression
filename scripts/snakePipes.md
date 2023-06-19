# snakePipes

The commands use a SNPfile & NMaskedIndex file. These can be created using the prepgenomes snakemake file in this repository (explained in the readme)

| dataset | SNPfile | NMASKINDEX_STAR | NMASKINDEX_BT | 
| --- | --- | --- | --- |
| 9sCa | all_CAST_EiJ_SNPs_129S1_SvImJ_reference.based_on_GRCm38.txt | STAR_9sCa/Genome | bt2_9sCa/genome |
| BlCa / CaBl / CaBlCaXi | all_SNPs_CAST_EiJ_GRCm38.txt.gz | STAR_BlCa/Genome | bt2_BlCa/genome |

Note that the organism argument (mm10_as) is created by createIndices previously (explained in the readme).

The sampleSheets used can be found under the samplesheets folder.
Other variables:

 - FQDIR: directory containing the RNAseq fastq files  
 - ODIR: output directory (will be created)  

## RNAseq runs

  > mRNA-seq -i $FQDIR -o $ODIR --libraryType 2 -j 20 --DAG --trim --trimmer trimgalore --trimmerOptions '-q 20 --fastqc --trim-n --clip_R1 6' -m allelic-mapping --SNPfile $SNPFILE --NMaskedIndex $NMASKINDEX_STAR --alignerOptions '--outFilterMultimapScoreRange 0 --limitBAMsortRAM 10000000000' --samplesheet $SAMPLESHEET mm10_as


  > Rscript scripts/allele_specific.degenes.cleanup.R $ODIR $FQDIR $SAMPLESHEET
  > Rscript scripts/allele_specific.degenes.category.R $ODIR $FQDIR $SAMPLESHEET

## TTseq

  > mRNA-seq -i $FQDIR -o $ODIR --libraryType 1 -m allelic-mapping --SNPfile $SNPFILE --NMaskedIndex $NMASKINDEX_STAR --sampleSheet $SAMPLESHEET mm10_as


## ChIPseq

The ChIPdicts are specified in the samplesheets folder and use either H3 or Input samples as controls.

  > DNA-mapping -i $FQDIR -o $ODIR -j 10 --DAG --fastqc --dedup --properPairs --mapq 3 --bwBinSize 25 --trim --trimmer trimgalore --trimmerOptions '-q 20 --fastqc --trim-n' -m allelic-mapping --SNPfile $SNPFILE --NMaskedIndex $NMASKINDEX_BT mm10_as

  > ChIP-seq -d $ODIR -j 10 mm10_as $CHIPDICT.yaml

Differential calling was performed with samplesheets as well, though for clarity only one example is giving in the repository:

 - ChIP_basic.tsv - bulk Control vs MSL2KO
 - ChIP_genome1.tsv - allele-specific Control vs MSL2KO
 - ChIP_genome2.tsv - allele-specific Control vs MSL2KO
 - ChIP_Control_g2vg1.tsv - genome 1 vs genome 2 for control samples
 - ChIP_MSL2KO1_g2vg1.tsv - genome 1 vs genome 2 for MSL2KO1 samples
 - ChIP_MSL2KO2_g2vg1.tsv - genome 1 vs genome 2 for MSL2KO2 samples

The calls slightly differ per mark:

For narrow marks:

  > ChIP-seq -d $ODIR -j 10 --windowSize 100 --FDR 1 --LFC 0 --CSAWCountsMethod Peak --sampleSheet $SAMPLESHEET mm10_as $CHIPDICT.yaml

For broad marks:

  > ChIP-seq -d $ODIR -j 10 --windowSize 500 --FDR 1 --LFC 0 --CSAWCountsMethod Peak --sampleSheet $SAMPLESHEET mm10_as $CHIPDICT.yaml

## ATACseq

  > DNA-mapping -i $FQDIR -o $ODIR -j 10 --DAG --fastqc --dedup --properPairs --mapq 3 --bwBinSize 25 --trim --trimmer trimgalore --trimmerOptions '-q 20 --fastqc --trim-n' -m allelic-mapping --SNPfile $SNPFILE --NMaskedIndex $NMASKINDEX_BT mm10_as

  > ATAC-seq -d $ODIR -j 10 --FDR 1 --LFC 0 --CSAWCountsMethod Peak --sampleSheet $SAMPLESHEET mm10