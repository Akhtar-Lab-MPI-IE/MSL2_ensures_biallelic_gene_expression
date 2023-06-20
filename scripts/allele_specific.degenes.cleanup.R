args = commandArgs(trailingOnly=TRUE)

prewd = args[1] #Output directory 
fastq.prewd = args[2] #Fastq directory
deseq2.directory = args[3] # SampleSheet
delete.chr = args[4]
chr.threshold = 1.2

library(rtracklayer)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tibble)
library(gridExtra)
library(ggpubr)


##Genome information from snakePipes pipeline
species.database.anno.gene = import(file.path(prewd,"Annotation/genes.filtered.gtf"))
exclude.category = c("IG_C_gene","IG_C_pseudogene","IG_D_gene","IG_D_pseudogene","IG_J_gene","IG_LV_gene","IG_pseudogene","IG_V_gene","IG_V_pseudogene","miRNA","misc_RNA","Mt_rRNA","Mt_tRNA","polymorphic_pseudogene","processed_pseudogene","pseudogene","ribozyme","rRNA","scaRNA","scRNA","snoRNA","snRNA","sRNA","TEC","TR_C_gene","TR_D_gene","TR_J_gene","TR_J_pseudogene","TR_V_gene","TR_V_pseudogene","transcribed_processed_pseudogene","transcribed_unitary_pseudogene","transcribed_unprocessed_pseudogene","translated_unprocessed_pseudogene","unitary_pseudogene","unprocessed_pseudogene") #c("snRNA", "processed_pseudogene", "lincRNA", "snoRNA", "misc_RNA", "rRNA", "scaRNA", "sRNA", "scRNA", "Mt_tRNA", "Mt_rRNA")
species.database.anno.gene.excluded = unique(species.database.anno.gene[species.database.anno.gene$gene_type %in% exclude.category,"gene_id"])
species.database.anno.gene.included = unique(species.database.anno.gene[species.database.anno.gene$gene_type %in% c("protein_coding"),"gene_id"])
mcols(species.database.anno.gene) = mcols(species.database.anno.gene)[,c("gene_id", "gene_name")]
id.name.correspondence = unique(as.data.frame(mcols(species.database.anno.gene)))
id.name.correspondence$gene_id = gsub("[.].*", "",id.name.correspondence$gene_id)
species.database.anno.gene.granges = unlist(reduce(split(species.database.anno.gene, species.database.anno.gene$gene_id), min.gapwidth=100000000))
names(species.database.anno.gene.granges) = sapply(strsplit(names(species.database.anno.gene.granges),"[.]"),function(x){x[1]})
species.database.anno.gene.granges$gene_id = names(species.database.anno.gene.granges)
species.database.anno.gene.granges$name = names(species.database.anno.gene.granges)
species.database.anno.gene.granges = species.database.anno.gene.granges[seqnames(species.database.anno.gene.granges) %in% c(1:19,"X","Y","MT")]
seqlevels(species.database.anno.gene.granges) =  c(1:19,"X","Y","MT")
species.database.anno.gene.granges = unique(species.database.anno.gene.granges)
lost.chr = GRanges(seqnames = 4, ranges = IRanges(start = 120000000, end = 156508116))
fdr = 0.01
lfcThreshold = 0
pvaluecol="pvalue" #"padj" #
basic.objects = c("basic.objects", ls(), "delete.chr")

if(TRUE){
  ##load in the original allelic DEseq file
  load(file.path(prewd, deseq2.directory,"DEseq_allelic_DESeq.Rdata"))
  dds.new=dds #[!rownames(dds) %in% species.database.anno.gene.granges.chrX]
  col.data.all = dds@colData
  col.data.all$combine = paste0(col.data.all$allele,col.data.all$condition) ##combine two conditions together
  colData(dds.new) = col.data.all
  dds.new = DESeqDataSet(dds.new, ~ combine ) ##change the design
  # dds.new <- estimateSizeFactors(dds.new)
  dds.new = DESeq(dds.new)
  print(resultsNames(dds.new))
  
  #compare genome1 and genome2 of control
  # dds.new=dds.new[!rownames(dds.new) %in% species.database.anno.gene.granges.chrX]
  # dds.new = estimateSizeFactors(dds.new)
  # dds.new = DESeq(dds.new)
  ddr.ck.alleles.comparison = results(dds.new, contrast = c("combine","genome2control","genome1control"), alpha = fdr,independentFiltering = TRUE)#foldchange between genome1 and genome2 
  ddr.ck.alleles.comparison =  as.data.frame(ddr.ck.alleles.comparison) #lfcShrink(dds.new, contrast = c("combine","genome2control","genome1control"), res = ddr.ck.alleles.comparison, type='normal')
  rownames(ddr.ck.alleles.comparison) = sapply(strsplit(rownames(ddr.ck.alleles.comparison),"[.]"),function(x){x[1]})
  ddr.ck.alleles.comparison$Status <- ifelse(is.na(ddr.ck.alleles.comparison[,pvaluecol]), "None", ifelse(ddr.ck.alleles.comparison[,pvaluecol] < fdr, ifelse(ddr.ck.alleles.comparison$log2FoldChange > 0, "UP", "DOWN"), "None"))
  ddr.ck.alleles.comparison.up = table(seqnames(species.database.anno.gene.granges[species.database.anno.gene.granges$gene_id %in% rownames(ddr.ck.alleles.comparison[ddr.ck.alleles.comparison$Status=="UP",])]))
  ddr.ck.alleles.comparison.down = table(seqnames(species.database.anno.gene.granges[species.database.anno.gene.granges$gene_id %in% rownames(ddr.ck.alleles.comparison[ddr.ck.alleles.comparison$Status=="DOWN",])]))
  print("genome2control / genome1control: up & down")
  print(ddr.ck.alleles.comparison.up)
  print(ddr.ck.alleles.comparison.down)
  
  ##plot the genome2 DEgenes number / genome1
  ratio = log2(ddr.ck.alleles.comparison.up / ddr.ck.alleles.comparison.down) 
  print("genome2control / genome1control up/down number log2ratio, should be less than 1")
  print(ratio)
  p = ggplot(as.data.frame(ratio), aes(x=Var1, y=Freq)) +
    geom_point() +
    ggtitle("CK") +
    geom_hline(yintercept = chr.threshold) + 
    geom_hline(yintercept = -chr.threshold) + 
    xlab("") + ylab("genome2 DEgenes number / genome1")+ theme_bw()
  
  #compare genome1 and genome2 of test
  ddr.ko.alleles.comparison = results(dds.new, contrast = c("combine","genome2test","genome1test"), alpha = fdr,independentFiltering = TRUE)#foldchange between genome1 and genome2 
  ddr.ko.alleles.comparison =  as.data.frame(ddr.ko.alleles.comparison) #lfcShrink(dds.new, contrast = c("combine","genome2test","genome1test"), res = ddr.ko.alleles.comparison, type='normal')
  rownames(ddr.ko.alleles.comparison) = sapply(strsplit(rownames(ddr.ko.alleles.comparison),"[.]"),function(x){x[1]})
  ddr.ko.alleles.comparison$Status <- ifelse(is.na(ddr.ko.alleles.comparison[,pvaluecol]), "None", ifelse(ddr.ko.alleles.comparison[,pvaluecol] < fdr, ifelse(ddr.ko.alleles.comparison$log2FoldChange > 0, "UP", "DOWN"), "None"))
  ddr.ko.alleles.comparison.up = table(seqnames(species.database.anno.gene.granges[species.database.anno.gene.granges$gene_id %in% rownames(ddr.ko.alleles.comparison[ddr.ko.alleles.comparison$Status=="UP",])]))
  ddr.ko.alleles.comparison.down = table(seqnames(species.database.anno.gene.granges[species.database.anno.gene.granges$gene_id %in% rownames(ddr.ko.alleles.comparison[ddr.ko.alleles.comparison$Status=="DOWN",])]))
  print("genome2test / genome1test: up & down")
  print(ddr.ko.alleles.comparison.up)
  print(ddr.ko.alleles.comparison.down)
  
  ##plot the genome2 DEgenes number / genome1
  ratio2 = log2(ddr.ko.alleles.comparison.up / ddr.ko.alleles.comparison.down) 
  print("genome2test / genome1test up/down number log2ratio, should be less than 1")
  print(ratio2)
  p2 = ggplot(as.data.frame(ratio2), aes(x=Var1, y=Freq)) +
    geom_point() +
    ggtitle("KO") +
    geom_hline(yintercept = chr.threshold) + 
    geom_hline(yintercept = -chr.threshold) + 
    xlab("") + ylab("genome2 DEgenes number / genome1")+ theme_bw()
  
  pdf(p2,file=file.path(prewd, deseq2.directory,"DEseq_resutls_CK_two_genomes_ratio_plot.png"))
  grid.arrange(p, p2)
  dev.off()
  
  ## determine the delete chromosomes based on the ratio < 1
  # delete.chr = c(3,12,11,14) 
  if(is.na(delete.chr)){
    delete.chr = as.numeric(setdiff( c(names(ratio[abs(ratio)>=chr.threshold]), names(ratio2[abs(ratio2)>=chr.threshold])) , c("X", "Y","MT") ))   #c(1,3,12,18)
    delete.chr = delete.chr[!is.na(delete.chr)]
  } else {
    delete.chr = unlist(strsplit(delete.chr, ","))
  }
  print("delete.chr:")
  print(delete.chr)
  species.database.anno.gene.granges.chrX = c() 
  if(length(delete.chr) != 0){
    if(any(delete.chr %in% "4")){species.database.anno.gene.granges.chrX = c(species.database.anno.gene.granges.chrX, unique(species.database.anno.gene.granges[queryHits(findOverlaps(species.database.anno.gene.granges,lost.chr))]$gene_id))}
    species.database.anno.gene.granges.chrX = c(species.database.anno.gene.granges.chrX, unique(species.database.anno.gene.granges[seqnames(species.database.anno.gene.granges) %in% setdiff(delete.chr,"4")]$gene_id))
    
  }
  print(sort(unique(seqnames(species.database.anno.gene.granges[species.database.anno.gene.granges.chrX]))))
  export(species.database.anno.gene.granges[species.database.anno.gene.granges.chrX], con = file.path(prewd, deseq2.directory,"delete.chr.bed"))
}

if(TRUE){
  ##deleted chromosomes
  species.database.anno.gene.granges.chrX = import(file.path(prewd, deseq2.directory,"delete.chr.bed"))$name

  ##load in the original allelic DEseq file
  load(file.path(prewd, deseq2.directory,"DEseq_allelic_DESeq.Rdata"))
  dds.new=dds[!gsub("[.].*", "", rownames(dds)) %in% species.database.anno.gene.granges.chrX]
  col.data.all = dds@colData
  col.data.all$combine = paste0(col.data.all$allele,col.data.all$condition) ##combine two conditions together
  colData(dds.new) = col.data.all
  dds.new = DESeqDataSet(dds.new, ~ combine ) ##change the design
  dds.new = DESeq(dds.new)
  print(resultsNames(dds.new))
  
  #compare genome1 and genome2 of control

  ddr.ck.alleles.comparison = results(dds.new, contrast = c("combine","genome2control","genome1control"), alpha = fdr,independentFiltering = TRUE)#foldchange between genome1 and genome2 
  ddr.ck.alleles.comparison =  as.data.frame(ddr.ck.alleles.comparison) #lfcShrink(dds.new, contrast = c("combine","genome2control","genome1control"), res = ddr.ck.alleles.comparison, type='normal')
  rownames(ddr.ck.alleles.comparison) = sapply(strsplit(rownames(ddr.ck.alleles.comparison),"[.]"),function(x){x[1]})
  ddr.ck.alleles.comparison[rownames(ddr.ck.alleles.comparison) %in% species.database.anno.gene.granges.chrX,"padj"]  = 1
  ddr.ck.alleles.comparison[rownames(ddr.ck.alleles.comparison) %in% species.database.anno.gene.granges.chrX,"pvalue"]  = 1
  ddr.ck.alleles.comparison$Status <- ifelse(is.na(ddr.ck.alleles.comparison[,pvaluecol]), "None", ifelse(ddr.ck.alleles.comparison[,pvaluecol] < fdr, ifelse(ddr.ck.alleles.comparison$log2FoldChange > 0, "UP", "DOWN"), "None"))
  ddr.ck.alleles.comparison[abs(ddr.ck.alleles.comparison$log2FoldChange) < lfcThreshold | is.na(ddr.ck.alleles.comparison$log2FoldChange),"Status"] <- "None"
  ddr.ck.alleles.comparison.up = table(seqnames(species.database.anno.gene.granges[species.database.anno.gene.granges$gene_id %in% rownames(ddr.ck.alleles.comparison[ddr.ck.alleles.comparison$Status=="UP",])]))
  ddr.ck.alleles.comparison.down = table(seqnames(species.database.anno.gene.granges[species.database.anno.gene.granges$gene_id %in% rownames(ddr.ck.alleles.comparison[ddr.ck.alleles.comparison$Status=="DOWN",])]))
  
  #compare genome1 and genome2 of test
  ddr.ko.alleles.comparison = results(dds.new, contrast = c("combine","genome2test","genome1test"), alpha = fdr,independentFiltering = TRUE)#foldchange between genome1 and genome2 
  ddr.ko.alleles.comparison =  as.data.frame(ddr.ko.alleles.comparison) #lfcShrink(dds.new, contrast = c("combine","genome2test","genome1test"), res = ddr.ko.alleles.comparison, type='normal')
  rownames(ddr.ko.alleles.comparison) = sapply(strsplit(rownames(ddr.ko.alleles.comparison),"[.]"),function(x){x[1]})
  ddr.ko.alleles.comparison[rownames(ddr.ko.alleles.comparison) %in% species.database.anno.gene.granges.chrX,"padj"]  = 1
  ddr.ko.alleles.comparison[rownames(ddr.ko.alleles.comparison) %in% species.database.anno.gene.granges.chrX,"pvalue"]  = 1
  ddr.ko.alleles.comparison$Status <- ifelse(is.na(ddr.ko.alleles.comparison[,pvaluecol]), "None", ifelse(ddr.ko.alleles.comparison[,pvaluecol] < fdr, ifelse(ddr.ko.alleles.comparison$log2FoldChange > 0, "UP", "DOWN"), "None"))
  ddr.ko.alleles.comparison[abs(ddr.ko.alleles.comparison$log2FoldChange) < lfcThreshold | is.na(ddr.ko.alleles.comparison$log2FoldChange),"Status"] <- "None"
  ddr.ko.alleles.comparison.up = table(seqnames(species.database.anno.gene.granges[species.database.anno.gene.granges$gene_id %in% rownames(ddr.ko.alleles.comparison[ddr.ko.alleles.comparison$Status=="UP",])]))
  ddr.ko.alleles.comparison.down = table(seqnames(species.database.anno.gene.granges[species.database.anno.gene.granges$gene_id %in% rownames(ddr.ko.alleles.comparison[ddr.ko.alleles.comparison$Status=="DOWN",])]))

  #compare genome1 of control and test
  ddr.genome1 = results(dds.new, contrast = c("combine","genome1test","genome1control"), alpha = fdr,independentFiltering = TRUE) #genome1 foldchange between CK and KO
  ddr.genome1 = as.data.frame(ddr.genome1) #lfcShrink(dds.new, contrast = c("combine","genome1test","genome1control"), res = ddr.genome1, type='normal')
  rownames(ddr.genome1) = sapply(strsplit(rownames(ddr.genome1),"[.]"),function(x){x[1]})
  ddr.genome1[rownames(ddr.genome1) %in% species.database.anno.gene.granges.chrX,"padj"]  = 1
  ddr.genome1[rownames(ddr.genome1) %in% species.database.anno.gene.granges.chrX,"pvalue"]  = 1
  ddr.genome1$Status <- ifelse(is.na(ddr.genome1[,pvaluecol]), "None", ifelse(ddr.genome1[,pvaluecol] < fdr, ifelse(ddr.genome1$log2FoldChange > 0, "UP", "DOWN"), "None"))
  ddr.genome1[abs(ddr.genome1$log2FoldChange) < lfcThreshold | is.na(ddr.genome1$log2FoldChange),"Status"] <- "None"
  print("genome1test / genome1control: up & down")
  print(table(seqnames(species.database.anno.gene.granges[species.database.anno.gene.granges$gene_id %in% rownames(ddr.genome1[ddr.genome1$Status=="UP",])])))
  print(table(seqnames(species.database.anno.gene.granges[species.database.anno.gene.granges$gene_id %in% rownames(ddr.genome1[ddr.genome1$Status=="DOWN",])])))
  
  #compare genome2 of control and test
  ddr.genome2 = results(dds.new, contrast = c("combine","genome2test","genome2control"), alpha = fdr,independentFiltering = TRUE) #genome2 foldchange between CK and KO
  ddr.genome2 =  as.data.frame(ddr.genome2) #lfcShrink(dds.new, contrast = c("combine","genome2test","genome2control"), res = ddr.genome2, type='normal')
  rownames(ddr.genome2) = sapply(strsplit(rownames(ddr.genome2),"[.]"),function(x){x[1]})
  ddr.genome2[rownames(ddr.genome2) %in% species.database.anno.gene.granges.chrX,"padj"]  = 1
  ddr.genome2[rownames(ddr.genome2) %in% species.database.anno.gene.granges.chrX,"pvalue"]  = 1
  ddr.genome2$Status <- ifelse(is.na(ddr.genome2[,pvaluecol]), "None", ifelse(ddr.genome2[,pvaluecol] < fdr, ifelse(ddr.genome2$log2FoldChange > 0, "UP", "DOWN"), "None"))
  ddr.genome2[abs(ddr.genome2$log2FoldChange) < lfcThreshold | is.na(ddr.genome2$log2FoldChange),"Status"] <- "None"
  print("genome2test / genome2control: up & down")
  print(table(seqnames(species.database.anno.gene.granges[species.database.anno.gene.granges$gene_id %in% rownames(ddr.genome2[ddr.genome2$Status=="UP",])])))
  print(table(seqnames(species.database.anno.gene.granges[species.database.anno.gene.granges$gene_id %in% rownames(ddr.genome2[ddr.genome2$Status=="DOWN",])])))
  
  
  save(ddr.genome1,ddr.genome2,ddr.ck.alleles.comparison,ddr.ko.alleles.comparison, file=file.path(prewd, deseq2.directory,"DEseq_allelic_genome_separated.RData"))
  write.table(ddr.ck.alleles.comparison, file=file.path(prewd, deseq2.directory,"DEseq_resutls_CK_two_genomes.tsv"),sep="\t",quote = FALSE,col.names = NA)
  write.table(ddr.ko.alleles.comparison, file=file.path(prewd, deseq2.directory,"DEseq_resutls_KO_two_genomes.tsv"),sep="\t",quote = FALSE,col.names = NA)
  write.table(ddr.genome1, file=file.path(prewd, deseq2.directory,"DEseq_resutls_genome1.tsv"),sep="\t",quote = FALSE,col.names = NA)
  write.table(ddr.genome2, file=file.path(prewd, deseq2.directory,"DEseq_resutls_genome2.tsv"),sep="\t",quote = FALSE,col.names = NA)
  
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}   

if(TRUE){
  species.database.anno.gene.granges.chrX = import(file.path(prewd, deseq2.directory,"delete.chr.bed"))$name
  genename.used = read.delim(file.path(prewd, deseq2.directory,"DEseq_basic_DEresults_LFCshrunk.tsv"),row.names = 1, stringsAsFactors = FALSE)
  rownames(genename.used) = gsub("[.].*", "", rownames(genename.used))
  ##Allelic DEanalysis using interaction
  load(file.path(prewd, deseq2.directory,"DEseq_allelic_DESeq.Rdata"))
  dds.allelic=dds[!gsub("[.].*", "", rownames(dds)) %in% species.database.anno.gene.granges.chrX]
  dds.allelic <- estimateSizeFactors(dds.allelic)
  print(head(1/sizeFactors(dds.allelic)))
  dds.allelic <- DESeq(dds.allelic,betaPrior = FALSE)
  print(resultsNames(dds.allelic))
  ddr = results(dds.allelic, lfcThreshold = 0, altHypothesis = "greaterAbs", name = "allelegenome2.conditiontest", alpha = 0.05, independentFiltering = TRUE)  #contrast = list(c("condition_test_vs_control", "allelegenome2.conditiontest"))
  ##remove duplicated or missing chromosomes
  peak.allelic = as.data.frame(ddr)
  rownames(peak.allelic) = sapply(strsplit(rownames(peak.allelic),"[.]"),function(x){x[1]})
  peak.allelic[rownames(peak.allelic) %in% species.database.anno.gene.granges.chrX,"padj"]  = 1
  peak.allelic[rownames(peak.allelic) %in% species.database.anno.gene.granges.chrX,"pvalue"]  = 1
  peak.allelic[rownames(peak.allelic) %in% species.database.anno.gene.granges.chrX,"Status"]  = "None"
  peak.allelic$Status <- ifelse(is.na(peak.allelic[,"pvalue"]), "None",ifelse(peak.allelic[,"pvalue"] < 0.05,ifelse(peak.allelic$log2FoldChange > 0, "UP", "DOWN"),"None"))
  peak.allelic[abs(peak.allelic$log2FoldChange) < 0.5 | is.na(peak.allelic$log2FoldChange),"Status"] <- "None"
  peak.allelic$external_gene_name  = genename.used[rownames(peak.allelic),"external_gene_name"]
  write.table(peak.allelic,file= file.path(prewd, deseq2.directory,"DEseq_allelic_DEresults.tsv"),quote = FALSE,sep = "\t")
  save(dds.allelic,ddr, file= file.path(prewd, deseq2.directory,"DEseq_allelic_DESeq.clean.Rdata"))

  load(file.path(prewd, deseq2.directory,"DEseq_basic_DESeq.Rdata"))
  dds.basic=dds[!gsub("[.].*", "", rownames(dds)) %in% species.database.anno.gene.granges.chrX]
  dds.basic <- estimateSizeFactors(dds.basic)
  dds.basic <- DESeq(dds.basic,betaPrior = FALSE)
  ddr = results(dds.basic, lfcThreshold = 0, altHypothesis = "greaterAbs", contrast = c("condition", "test", "control"), alpha = fdr, independentFiltering = TRUE)  
  ddr = lfcShrink(dds.basic, contrast = c("condition", "test", "control"), res = ddr, type='normal')
  ##remove duplicated or missing chromosomes
  peak.basic = as.data.frame(ddr)
  rownames(peak.basic) = sapply(strsplit(rownames(peak.basic),"[.]"),function(x){x[1]})
  peak.basic[rownames(peak.basic) %in% species.database.anno.gene.granges.chrX,"padj"]  = 1
  peak.basic[rownames(peak.basic) %in% species.database.anno.gene.granges.chrX,"pvalue"]  = 1
  peak.basic[rownames(peak.basic) %in% species.database.anno.gene.granges.chrX,"Status"]  = "None"
  peak.basic$Status <- ifelse(is.na(peak.basic[,"padj"]), "None", ifelse(peak.basic[,"padj"] < fdr,ifelse(peak.basic$log2FoldChange > 0, "UP", "DOWN"),"None"))
  peak.basic[abs(peak.basic$log2FoldChange) < lfcThreshold | is.na(peak.basic$log2FoldChange),"Status"] <- "None"
  peak.basic$external_gene_name  = genename.used[rownames(peak.basic),"external_gene_name"]
  write.table(peak.basic,file= file.path(prewd, deseq2.directory,"DEseq_basic_DEresults.tsv"),quote = FALSE,sep = "\t")
  save(dds.basic, ddr, file= file.path(prewd, deseq2.directory,"DEseq_basic_DESeq.clean.Rdata"))

  ###scale factor
  normal.normalization.factor = as.data.frame(1/sizeFactors(dds.basic))
  rownames(normal.normalization.factor) =  gsub("_all", "", rownames(normal.normalization.factor))
  write.table(normal.normalization.factor,file=file.path(fastq.prewd, "readme/normal_normalization.txt"),col.names = FALSE, quote = FALSE,sep=":")
  allelic.normalization.factor = as.data.frame(1/sizeFactors(dds.allelic))
  rownames(allelic.normalization.factor) = gsub("_genome", ".genome", rownames(allelic.normalization.factor))
  write.table(allelic.normalization.factor,file=file.path(fastq.prewd, "readme/allelic_normalization.txt"),col.names = FALSE, quote = FALSE,sep=":")

  rm(list=setdiff(ls(),basic.objects))
  gc()
}
