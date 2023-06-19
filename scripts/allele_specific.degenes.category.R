##Setup#####
args = commandArgs(trailingOnly=TRUE)
prewd = args[1] #Output directory
output.dir= args[2] #Fastq directory
deseq2.directory = args[3] #SampleSheet

library(rtracklayer)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tibble)
library(gridExtra)
library(ggpubr)
library(pheatmap)

##Genome information from snakePipes pipeline
if(TRUE){
  species.database.anno.gene = import(file.path(prewd,"Annotation/genes.filtered.gtf")) #import("/data/repository/organisms/GRCm38_ensembl/Ensembl/release-91/genes.gtf") #
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
  gene.width = as.data.frame(width(species.database.anno.gene.granges) ) 
  rownames(gene.width) = species.database.anno.gene.granges$gene_id
  species.database.anno.gene.granges = species.database.anno.gene.granges[seqnames(species.database.anno.gene.granges) %in% c(1:19,"X","Y","MT")]
  seqlevels(species.database.anno.gene.granges) =  c(1:19,"X","Y","MT")
  species.database.anno.gene.granges = unique(species.database.anno.gene.granges)
}

fdr = 0.01 ##The outside cutoff, should be same as cleanup script
pvaluecol="pvalue" #"padj" #"pvalue" #
basic.objects = c("basic.objects",ls())
dir.create(file.path(output.dir))
dir.create(file.path(output.dir,"Expression"))

##Allelic RNA-seq expression level UP and Down bed file
if(TRUE){
  
  ##allelic DE genes
  peak.allelic = read.delim(file.path(prewd, deseq2.directory,"DEseq_allelic_DEresults.tsv"),row.names = 1)
  peak.allelic.up = unlist(sapply(strsplit(rownames(peak.allelic[peak.allelic$Status=="UP",]),"[.]"),function(x){x[1]})) 
  peak.allelic.down = unlist(sapply(strsplit(rownames(peak.allelic[peak.allelic$Status=="DOWN",]),"[.]"),function(x){x[1]})) 
  print(paste0("allelic up genes are ",length(peak.allelic.up)))
  print(paste0("allelic down genes are ",length(peak.allelic.down)))
  
  peak.allelic.up.genes = species.database.anno.gene.granges[species.database.anno.gene.granges$gene_id %in% peak.allelic.up]
  peak.allelic.down.genes = species.database.anno.gene.granges[species.database.anno.gene.granges$gene_id %in% peak.allelic.down]
  print(paste0("check: export up genes bed file are ",length(peak.allelic.up.genes)))
  print(paste0("check: export down genes bed file are ",length(peak.allelic.down.genes)))
  
  export(granges(peak.allelic.up.genes),con=file.path(prewd, deseq2.directory,"DE_allelic_significant_up_genes.bed")) #[!as.character(seqnames(peak.allelic.up.genes)) %in% c("3","12")]
  export(granges(peak.allelic.down.genes),con=file.path(prewd, deseq2.directory,"DE_allelic_significant_down_genes.bed")) #[!as.character(seqnames(peak.allelic.down.genes)) %in% c("3","12")]
  
  ##basic DE genes
  peak.basic = read.delim(file.path(prewd, deseq2.directory,"DEseq_basic_DEresults.tsv"),row.names = 1)
  peak.basic.up = unlist(sapply(strsplit(rownames(peak.basic[peak.basic$Status=="UP",]),"[.]"),function(x){x[1]}))
  peak.basic.down = unlist(sapply(strsplit(rownames(peak.basic[peak.basic$Status=="DOWN",]),"[.]"),function(x){x[1]}))
  print(paste0("basic up genes are ",length(peak.basic.up)))
  print(paste0("basic down genes are ",length(peak.basic.down)))
  
  peak.basic.up.genes = species.database.anno.gene.granges[species.database.anno.gene.granges$gene_id %in% peak.basic.up]
  peak.basic.down.genes = species.database.anno.gene.granges[species.database.anno.gene.granges$gene_id %in% peak.basic.down]
  print(paste0("check: export basic up genes bed file are ",length(peak.basic.up.genes)))
  print(paste0("check: export basic down genes bed file are ",length(peak.basic.down.genes)))
  
  export(granges(peak.basic.up.genes),con=file.path(prewd, deseq2.directory,"DE_basic_significant_up_genes.bed")) #[!as.character(seqnames(peak.basic.up.genes)) %in% c("3","12")]
  export(granges(peak.basic.down.genes),con=file.path(prewd, deseq2.directory,"DE_basic_significant_down_genes.bed")) #[!as.character(seqnames(peak.basic.down.genes)) %in% c("3","12")]
  
  #Select Random genes as negative control (exluding  both basic and allelic down and up genes)
  up.down.genes = unique(c(peak.allelic.up.genes$gene_id,peak.allelic.down.genes$gene_id,peak.basic.up.genes$gene_id,peak.basic.down.genes$gene_id))
  rownames(peak.allelic) = gsub("[.][0-9]*","",rownames(peak.allelic))
  peak.allelic.expressed = peak.allelic[(!(rownames(peak.allelic) %in% up.down.genes)),]
  peak.allelic.expressed = peak.allelic.expressed[peak.allelic.expressed$baseMean>2,]
  set.seed(100)
  peak.allelic.expressed.selected = peak.allelic.expressed[sample(nrow(peak.allelic.expressed),300),]
  export(species.database.anno.gene.granges[names(species.database.anno.gene.granges) %in% rownames(peak.allelic.expressed.selected)],con=file.path(output.dir,"random.genes.allelic.bed"))
  
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}

if(TRUE){
  ##foldchange allelic calculated by DESeq2###################################
  ########################################################################
  ##README: log2FC < 0 means genome1 is expressed while genome2 down regulated
  ddr.allelic = read.delim(file.path(prewd, deseq2.directory,"DEseq_allelic_DEresults.tsv"))
  rownames(ddr.allelic) = sapply(strsplit(rownames(ddr.allelic),"[.]"),function(x){x[1]})
  ddr.basic = read.delim(file.path(prewd, deseq2.directory,"DEseq_basic_DEresults.tsv"))
  rownames(ddr.basic) = sapply(strsplit(rownames(ddr.basic),"[.]"),function(x){x[1]})
  colnames(ddr.basic) = paste0("basic_", colnames(ddr.basic))
  
  ###Count Info##########################################
  ########################################################
  #load basic counts information
  load(file.path(prewd, deseq2.directory,"DEseq_basic_DESeq.Rdata"))
  counts.basic = counts(dds,normalized=TRUE)
  rownames(counts.basic) = sapply(strsplit(rownames(counts.basic),"[.]"),function(x){x[1]})
  if(length(grep("^L_",colnames(counts.basic))) != 0){
    samplenames = unique(sapply(strsplit(colnames(counts.basic),"_"),function(x){x[7]}))
  }else{
    samplenames = unique(sapply(strsplit(colnames(counts.basic),"_"),function(x){x[2]}))
  }
  counts.basic.mean = cbind()
  for (samplename in samplenames){
    counts.basic.mean = cbind(counts.basic.mean,apply(counts.basic[,grep(samplename,colnames(counts.basic)),drop=FALSE],1,mean))
  }
  colnames(counts.basic.mean) = samplenames
  counts.basic.mean = as.data.frame(counts.basic.mean)
  # ####counts.basic.mean = counts.basic.mean * 1000 / gene.width[rownames(counts.basic.mean),1]
  
  ##load allelic counts information
  load(file.path(prewd, deseq2.directory,"DEseq_allelic_DESeq.Rdata"))
  counts.allelic = counts(dds,normalized=TRUE)
  rownames(counts.allelic) = sapply(strsplit(rownames(counts.allelic),"[.]"),function(x){x[1]})
  
  if(length(grep("^L_",colnames(counts.allelic))) != 0){
    colnames(counts.allelic) = sapply(strsplit(colnames(counts.allelic),"_"),function(x){paste(x[6],x[7],x[8],x[3],sep="_")})
    samplenames = unique(sapply(strsplit(colnames(counts.allelic),"_"),function(x){paste(x[2],x[3],sep="_")}))
  }else{
    colnames(counts.allelic) = sapply(strsplit(colnames(counts.allelic),"_"),function(x){paste(x[1],x[2],x[4],x[3],sep="_")})
    samplenames = unique(sapply(strsplit(colnames(counts.allelic),"_"),function(x){paste(x[2],x[3],sep="_")}))
  }
  counts.allelic.mean = cbind()
  for (samplename in samplenames){
    counts.allelic.mean = cbind(counts.allelic.mean,apply(counts.allelic[,grep(samplename,colnames(counts.allelic)),drop=FALSE],1,mean))
  }
  counts.allelic.mean = as.data.frame(counts.allelic.mean)
  # ####counts.allelic.mean = counts.allelic.mean * 1000 / gene.width[rownames(counts.allelic.mean),1]
  
  if(dim(counts.allelic.mean)[2] == 6){
    colnames(counts.basic.mean) = c("CK","KO1","KO2")
    colnames(counts.allelic.mean) = c("CK_genome1","CK_genome2","KO1_genome1","KO1_genome2","KO2_genome1","KO2_genome2")
    counts.allelic.mean$CK_genome1_ratio = counts.allelic.mean$CK_genome1 / (counts.allelic.mean$CK_genome1 + counts.allelic.mean$CK_genome2)
    counts.allelic.mean$KO1_genome1_ratio = counts.allelic.mean$KO1_genome1 / (counts.allelic.mean$KO1_genome1 + counts.allelic.mean$KO1_genome2)
    counts.allelic.mean$KO2_genome1_ratio = counts.allelic.mean$KO2_genome1 / (counts.allelic.mean$KO2_genome1 + counts.allelic.mean$KO2_genome2)
  }  else if (dim(counts.allelic.mean)[2] == 4){
    colnames(counts.basic.mean) = c("CK","KO")
    colnames(counts.allelic.mean) = c("CK_genome1","CK_genome2","KO_genome1","KO_genome2")
    counts.allelic.mean$CK_genome1_ratio = counts.allelic.mean$CK_genome1 / (counts.allelic.mean$CK_genome1 + counts.allelic.mean$CK_genome2)
    counts.allelic.mean$KO1_genome1_ratio = counts.allelic.mean$KO_genome1 / (counts.allelic.mean$KO_genome1 + counts.allelic.mean$KO_genome2)
   }
  
  counts.normalized.separated.genome.allgenes = cbind(counts.basic.mean[rownames(counts.allelic.mean),], counts.allelic.mean)
  write.table(counts.normalized.separated.genome.allgenes,file = file.path(output.dir,"Expression","counts.normalized.separated.genome.allgenes.txt"), sep="\t",quote = FALSE, col.names = NA)
  save(counts.normalized.separated.genome.allgenes,counts.basic.mean,counts.allelic.mean,file = file.path(output.dir,"Expression","counts.normalized.separated.genome.allgenes.RData"))
  
  ###foldchange calculation for genome1 and genome2 seperately using DEseq2########################
  ###foldchange are shrinked#########################################################
  load(file.path(prewd, deseq2.directory,"DEseq_allelic_genome_separated.RData"))
  
  ###calculation of log2FC seperating the genome by my own formula ##########
  ########################################################################
  ##separate foldchange for genome1 and genome2, calculated by myself
  ##README: If WT and KO are both 0, the result will be NaN; if only KO is 0, it will be -Inf; if only WT is 0, it will be Inf; KO1 foldchange = log2( (KO1_G2/WT_G2)/(KO1_G1/WT_G1) ): -1 (-Inf) means G1_biased; 1(Inf)
  if(dim(counts.allelic.mean)[2] == 9){
    counts.allelic.mean.log2FC = as.data.frame(counts.allelic.mean + 1)
    colnames(counts.allelic.mean.log2FC) = paste(colnames(counts.allelic.mean.log2FC), "log2FC",sep="_")
    #my own calculation
    counts.allelic.mean.log2FC$KO1_log2FC = log2((counts.allelic.mean.log2FC[,"KO1_genome2_log2FC"]/counts.allelic.mean.log2FC[,"CK_genome2_log2FC"])/(counts.allelic.mean.log2FC[,"KO1_genome1_log2FC"]/counts.allelic.mean.log2FC[,"CK_genome1_log2FC"]))
    counts.allelic.mean.log2FC$KO2_log2FC = log2((counts.allelic.mean.log2FC[,"KO2_genome2_log2FC"]/counts.allelic.mean.log2FC[,"CK_genome2_log2FC"])/(counts.allelic.mean.log2FC[,"KO2_genome1_log2FC"]/counts.allelic.mean.log2FC[,"CK_genome1_log2FC"]))
    counts.allelic.mean.log2FC$Mean_log2FC = apply(counts.allelic.mean.log2FC[,c("KO1_log2FC","KO2_log2FC")],1,mean)
    counts.allelic.mean.log2FC$Std_log2FC = apply(counts.allelic.mean.log2FC[,c("KO1_log2FC","KO2_log2FC")],1,function(x) sd(x)/sqrt(length(x)))
    counts.allelic.mean.log2FC[,grep("genome1",colnames(counts.allelic.mean.log2FC))] = log2(counts.allelic.mean.log2FC[,grep("genome1",colnames(counts.allelic.mean.log2FC))]/counts.allelic.mean.log2FC[,"CK_genome1_log2FC"])
    counts.allelic.mean.log2FC[,grep("genome2",colnames(counts.allelic.mean.log2FC))] = log2(counts.allelic.mean.log2FC[,grep("genome2",colnames(counts.allelic.mean.log2FC))]/counts.allelic.mean.log2FC[,"CK_genome2_log2FC"])
    counts.allelic.mean.log2FC = counts.allelic.mean.log2FC[,-c(1:2)]
    #DEseq2 calculation
    counts.allelic.mean.log2FC$KOAll_genome1_log2FC = ddr.genome1[rownames(counts.allelic.mean.log2FC),"log2FoldChange"]
    counts.allelic.mean.log2FC$KOAll_genome1_pvalue = ddr.genome1[rownames(counts.allelic.mean.log2FC),pvaluecol]
    counts.allelic.mean.log2FC$KOAll_genome2_log2FC = ddr.genome2[rownames(counts.allelic.mean.log2FC),"log2FoldChange"]
    counts.allelic.mean.log2FC$KOAll_genome2_pvalue = ddr.genome2[rownames(counts.allelic.mean.log2FC),pvaluecol]
    counts.allelic.mean.log2FC$AllelicComparison_log2FC = ddr.ck.alleles.comparison[rownames(counts.allelic.mean.log2FC),"log2FoldChange"]
    counts.allelic.mean.log2FC$AllelicComparison_pvalue = ddr.ck.alleles.comparison[rownames(counts.allelic.mean.log2FC),pvaluecol]
    counts.allelic.mean.log2FC$AllelicComparisonKO_log2FC = ddr.ko.alleles.comparison[rownames(counts.allelic.mean.log2FC),"log2FoldChange"]
    counts.allelic.mean.log2FC$AllelicComparisonKO_pvalue = ddr.ko.alleles.comparison[rownames(counts.allelic.mean.log2FC),pvaluecol]
    colnames.order = c(grep("genome1",colnames(counts.allelic.mean.log2FC)),grep("genome2",colnames(counts.allelic.mean.log2FC)),grep("genome",colnames(counts.allelic.mean.log2FC),invert = TRUE))
    counts.allelic.mean.log2FC = counts.allelic.mean.log2FC[,colnames.order]
  } else if (dim(counts.allelic.mean)[2] == 6){
    counts.allelic.mean.log2FC = as.data.frame(counts.allelic.mean + 1)
    colnames(counts.allelic.mean.log2FC) = paste(colnames(counts.allelic.mean.log2FC), "log2FC",sep="_")
    #my own calculation
    counts.allelic.mean.log2FC$KO1_log2FC = log2((counts.allelic.mean.log2FC[,"KO_genome2_log2FC"]/counts.allelic.mean.log2FC[,"CK_genome2_log2FC"])/(counts.allelic.mean.log2FC[,"KO_genome1_log2FC"]/counts.allelic.mean.log2FC[,"CK_genome1_log2FC"]))
    counts.allelic.mean.log2FC[,grep("genome1",colnames(counts.allelic.mean.log2FC))] = log2(counts.allelic.mean.log2FC[,grep("genome1",colnames(counts.allelic.mean.log2FC))]/counts.allelic.mean.log2FC[,"CK_genome1_log2FC"])
    counts.allelic.mean.log2FC[,grep("genome2",colnames(counts.allelic.mean.log2FC))] = log2(counts.allelic.mean.log2FC[,grep("genome2",colnames(counts.allelic.mean.log2FC))]/counts.allelic.mean.log2FC[,"CK_genome2_log2FC"])
    counts.allelic.mean.log2FC = counts.allelic.mean.log2FC[,-c(1:2)]
    #DEseq2 calculation
    counts.allelic.mean.log2FC$KOAll_genome1_log2FC = ddr.genome1[rownames(counts.allelic.mean.log2FC),"log2FoldChange"]
    counts.allelic.mean.log2FC$KOAll_genome1_pvalue = ddr.genome1[rownames(counts.allelic.mean.log2FC),pvaluecol]
    counts.allelic.mean.log2FC$KOAll_genome2_log2FC = ddr.genome2[rownames(counts.allelic.mean.log2FC),"log2FoldChange"]
    counts.allelic.mean.log2FC$KOAll_genome2_pvalue = ddr.genome2[rownames(counts.allelic.mean.log2FC),pvaluecol]
    counts.allelic.mean.log2FC$AllelicComparison_log2FC = ddr.ck.alleles.comparison[rownames(counts.allelic.mean.log2FC),"log2FoldChange"]
    counts.allelic.mean.log2FC$AllelicComparison_pvalue = ddr.ck.alleles.comparison[rownames(counts.allelic.mean.log2FC),pvaluecol]
    counts.allelic.mean.log2FC$AllelicComparisonKO_log2FC = ddr.ko.alleles.comparison[rownames(counts.allelic.mean.log2FC),"log2FoldChange"]
    counts.allelic.mean.log2FC$AllelicComparisonKO_pvalue = ddr.ko.alleles.comparison[rownames(counts.allelic.mean.log2FC),pvaluecol]
  }
  
  ##combine my calculation for log2FC and DEseq2 log2FC, mean counts for each sample together
  ########################################################################################################
  counts.allelic.mean.log2FC.combined = cbind(counts.allelic.mean[rownames(ddr.allelic),], counts.allelic.mean.log2FC[rownames(ddr.allelic),], ddr.basic[rownames(ddr.allelic),c(1,2,6,7)], ddr.allelic[,-c(3:4)])
  counts.allelic.mean.log2FC.combined[is.na(counts.allelic.mean.log2FC.combined$log2FoldChange), "log2FoldChange"] = 0
  counts.allelic.mean.log2FC.combined[is.na(counts.allelic.mean.log2FC.combined[,"pvalue"]), "pvalue"] = 1
  counts.allelic.mean.log2FC.combined[is.na(counts.allelic.mean.log2FC.combined[,"padj"]), "padj"] = 1
  counts.allelic.mean.log2FC.combined[is.na(counts.allelic.mean.log2FC.combined$KOAll_genome1_log2FC), "KOAll_genome1_log2FC"] = 0
  counts.allelic.mean.log2FC.combined[is.na(counts.allelic.mean.log2FC.combined$KOAll_genome1_pvalue), "KOAll_genome1_pvalue"] = 1
  counts.allelic.mean.log2FC.combined[is.na(counts.allelic.mean.log2FC.combined$KOAll_genome2_log2FC), "KOAll_genome2_log2FC"] = 0
  counts.allelic.mean.log2FC.combined[is.na(counts.allelic.mean.log2FC.combined$KOAll_genome2_pvalue), "KOAll_genome2_pvalue"] = 1
  counts.allelic.mean.log2FC.combined$Quantile = "log2FC > -1"
  write.table(counts.allelic.mean.log2FC.combined,file = file.path(output.dir,"Expression","log2FC.separated.genome.allgenes.txt"), sep="\t",quote = FALSE, col.names = NA)
  
  log2FC.separated.genome.allgenes = counts.allelic.mean.log2FC.combined
  save(log2FC.separated.genome.allgenes,file = file.path(output.dir,"Expression","log2FC.separated.genome.allgenes.RData"))
  
  # rm(list=setdiff(ls(),basic.objects))
  # gc()
}

if(TRUE){
  ##extract DEgenes identified in both basic and allelic analysis
  ########################################################################
  ##load basic and allelic DEgenes info
  allelic.degenes = read.delim(file.path(prewd, deseq2.directory,"DEseq_allelic_DEresults.tsv"), row.names = 1, stringsAsFactors = FALSE)
  allelic.degenes[is.na(allelic.degenes$log2FoldChange), "log2FoldChange"] = 0
  allelic.degenes[is.na(allelic.degenes[,pvaluecol]), pvaluecol] = 1
  basic.degenes = read.delim(file.path(prewd, deseq2.directory,"DEseq_basic_DEresults.tsv"), row.names = 1, stringsAsFactors = FALSE)
  basic.degenes[is.na(basic.degenes$log2FoldChange), "log2FoldChange"] = 0
  basic.degenes[is.na(basic.degenes[,pvaluecol]), pvaluecol] = 1
  # combined.degenesfdr = unique(c(rownames(basic.degenes[basic.degenes$log2FoldChange <0 & basic.degenes[,pvaluecol]<=0.01,]), rownames(allelic.degenes[allelic.degenes[,pvaluecol]<=fdr,])) )
  combined.degenes0.1 = unique(c(rownames(basic.degenes[basic.degenes$log2FoldChange < 0 & basic.degenes[,"padj"]<=fdr,]), rownames(allelic.degenes[allelic.degenes[,pvaluecol] <= 0.05 & abs(allelic.degenes$log2FoldChange) >= 0.5,])) )
  # combined.degenes0.1.new = unique(c(import(file.path(prewd, deseq2.directory,"DE_basic_significant_down_genes.bed")), import(file.path(prewd, deseq2.directory,"DE_allelic_significant_down_genes.bed")), import(file.path(prewd, deseq2.directory,"DE_allelic_significant_up_genes.bed")))$name)
  ##remove upregulated genes and non-allelic genes
  load(file.path(output.dir,"Expression","log2FC.separated.genome.allgenes.RData"))
  counts.allelic.mean.log2FC.combined = log2FC.separated.genome.allgenes[!rownames(log2FC.separated.genome.allgenes) %in% gsub("[.].*", "", species.database.anno.gene.excluded$gene_id),]
  # counts.allelic.mean.log2FC.combined.degenes = counts.allelic.mean.log2FC.combined[rownames(counts.allelic.mean.log2FC.combined) %in% combined.degenes0.1,]
  counts.allelic.mean.log2FC.combined.degenes = counts.allelic.mean.log2FC.combined %>% filter( (KOAll_genome1_log2FC < 0 & KOAll_genome1_pvalue <= fdr) | (KOAll_genome2_log2FC < 0 & KOAll_genome2_pvalue <= fdr) | rownames(counts.allelic.mean.log2FC.combined) %in% combined.degenes0.1 )
  counts.allelic.mean.log2FC.combined.degenes = counts.allelic.mean.log2FC.combined.degenes[order(counts.allelic.mean.log2FC.combined.degenes$log2FoldChange),]
  print("all basic and allelic DEgene number")
  print(nrow(counts.allelic.mean.log2FC.combined.degenes))
  print("remove upregulated genes (G1,G2 log2FC all up )")
  tmp = counts.allelic.mean.log2FC.combined.degenes[,c(grep("genome1_log2FC",colnames(counts.allelic.mean.log2FC.combined.degenes)),grep("genome2_log2FC",colnames(counts.allelic.mean.log2FC.combined.degenes)))]
  counts.allelic.mean.log2FC.combined.degenes = counts.allelic.mean.log2FC.combined.degenes[!apply(tmp,1,function(x){return(all(x>0))}),] 
  # #remove basic up regulated genes, this step is already deleted
  # basic.up.genes = import(file.path(prewd, deseq2.directory,"DE_basic_significant_up_genes.bed"))
  # counts.allelic.mean.log2FC.combined.degenes = counts.allelic.mean.log2FC.combined.degenes[!rownames(counts.allelic.mean.log2FC.combined.degenes) %in% basic.up.genes$name,]
  print(nrow(counts.allelic.mean.log2FC.combined.degenes))
  print("remove non-expressed genes (basic expressed but got no SNPs in allele specific)")
  tmp = counts.allelic.mean.log2FC.combined.degenes[,1:4]
  tmp = tmp[,grep("_log2FC",colnames(tmp),invert = TRUE)]
  counts.allelic.mean.log2FC.combined.degenes = counts.allelic.mean.log2FC.combined.degenes[!apply(tmp,1,function(x){return(all(x==0))}),] 
  print(nrow(counts.allelic.mean.log2FC.combined.degenes))
  
  ########category of all DEgenes ##########all:either allelic DE analysis pvalue <= 0.1 and g2 down or single allele DE analysis pvalue <= fdr and g2 down
  ##########################################normal: logfc0.5, g2 down > 0.5; logfc1: logfc1, g2 down > 1; logfc1.extend: logfc1, g2 down > 1 and allelic DE analysis pvalue <= 0.1
  ##########################################extend:allelic DE analysis pvalue <= 0.1; dim(one allele down but the other one up)
  counts.allelic.mean.log2FC.combined.degenes = counts.allelic.mean.log2FC.combined.degenes[order(abs(counts.allelic.mean.log2FC.combined.degenes$log2FoldChange),decreasing = TRUE),]
  
  pvalue.cutoff= 0.05  ##Inner cutoff for allelic level
  pvalue.cutoff.loose= 0.05  ##Inner cutoff for allelic level
  pvalue.cutoff.bi= as.numeric(args[4])  ##also 0.05
  log2foldchange.lowest=-0.2 #-0.5 ##allelic lfc threshold
  log2foldchange.small=1 #0.5 ##allelic lfc threshold
  log2foldchange=2 #1 ##allelic lfc threshold
  log2foldchange.bi=log2foldchange #as.numeric(args[5]) ##allelic lfc threshold
  log2foldchange.big = 3 #2 # 
  mono.log2fc = 1 #2.5 #  mono.log2fc>0 means G2 expressed more, and mono.log2fc < 0 means G1 expressed more
  # logfc1.fdr= 0.1  ##fdr cutoff only for logfc >=1 degenes
  # expression.cutoff = 0.1
  ###########
  ##Bi to mono
  ###########
  counts.allelic.mean.log2FC.degenes.bi = counts.allelic.mean.log2FC.combined.degenes #counts.allelic.mean.log2FC.combined.degenes[!rownames(counts.allelic.mean.log2FC.combined.degenes) %in% rownames(counts.allelic.mean.log2FC.degenes.mono.order),]
  
  ##G1 retained
  counts.allelic.mean.log2FC.degenes.bi.genome1.biased.extend = counts.allelic.mean.log2FC.degenes.bi %>% filter( (log2FoldChange < -0.5 & get(pvaluecol) <= pvalue.cutoff.loose) )  %>% filter( AllelicComparison_log2FC <= (mono.log2fc) )  %>% filter( KOAll_genome2_log2FC <= 0) #CK_genome1_ratio <= 0.9 & counts.allelic.mean.log2FC.degenes.bi[counts.allelic.mean.log2FC.degenes.bi$log2FoldChange < 0 & counts.allelic.mean.log2FC.degenes.bi[,pvaluecol] <= pvalue.cutoff.extend  & counts.allelic.mean.log2FC.degenes.bi$AllelicComparison_log2FC <= (mono.log2fc) ,] # & counts.allelic.mean.log2FC.degenes.bi$AllelicComparison_pvalue >= pvalue.cutoff #[counts.allelic.mean.log2FC.degenes.bi.genome1.biased$log2FoldChange <= -(log2foldchange),] #allelic foldchange <= -1
  counts.allelic.mean.log2FC.degenes.bi.genome1.biased = counts.allelic.mean.log2FC.degenes.bi %>% filter( (log2FoldChange < -0.5 & get(pvaluecol) <= pvalue.cutoff.bi) )  %>% filter( AllelicComparison_log2FC <= (mono.log2fc) )  %>% filter(CK_genome1_ratio <= 0.9 &  KOAll_genome2_log2FC <= log2foldchange.lowest) # & AllelicComparison_pvalue >= pvalue.cutoff Status=="DOWN", exclude G2 expressed more than 2 fold than G1, must #(KOAll_genome2_pvalue <= fdr & KOAll_genome2_log2FC < 0 & KOAll_genome1_log2FC * 5 > KOAll_genome2_log2FC & KOAll_genome1_log2FC > -log2foldchange.small) | 
  counts.allelic.mean.log2FC.degenes.bi.genome1.biased.logfc0.5 = counts.allelic.mean.log2FC.degenes.bi.genome1.biased  %>% filter( KOAll_genome2_log2FC < -(log2foldchange.small) )  #[apply(counts.allelic.mean.log2FC.degenes.bi.genome1.biased.logfc1[,grep("KOAll_genome2_log2FC",colnames(counts.allelic.mean.log2FC.degenes.bi.genome1.biased.logfc1)),drop=FALSE],1,function(x){any(x<= -(log2foldchange.small))}),] #G2 is down regulated 2 fold, must
  counts.allelic.mean.log2FC.degenes.bi.genome1.biased.logfc1 = counts.allelic.mean.log2FC.degenes.bi.genome1.biased  %>% filter( KOAll_genome2_log2FC < -(log2foldchange.bi) ) #[apply(counts.allelic.mean.log2FC.degenes.bi.genome1.biased.logfc1[,grep("KOAll_genome2_log2FC",colnames(counts.allelic.mean.log2FC.degenes.bi.genome1.biased.logfc1)),drop=FALSE],1,function(x){any(x<= -(log2foldchange))}),] #G2 is down regulated 2 fold, must
  counts.allelic.mean.log2FC.degenes.bi.genome1.biased.quantile.first = counts.allelic.mean.log2FC.degenes.bi.genome1.biased  %>% filter( KOAll_genome2_log2FC > -(log2foldchange.small) & KOAll_genome2_log2FC < 0 )
  counts.allelic.mean.log2FC.degenes.bi.genome1.biased.quantile.second = counts.allelic.mean.log2FC.degenes.bi.genome1.biased  %>% filter( KOAll_genome2_log2FC > -(log2foldchange.bi) & KOAll_genome2_log2FC <= -(log2foldchange.small) )
  counts.allelic.mean.log2FC.degenes.bi.genome1.biased.quantile.third = counts.allelic.mean.log2FC.degenes.bi.genome1.biased  %>% filter( KOAll_genome2_log2FC < -(log2foldchange.bi) )
  counts.allelic.mean.log2FC.degenes.bi.genome1.biased.quantile.thirdpart = counts.allelic.mean.log2FC.degenes.bi.genome1.biased  %>% filter( KOAll_genome2_log2FC > -(log2foldchange.big) & KOAll_genome2_log2FC <= -(log2foldchange.bi) )
  counts.allelic.mean.log2FC.degenes.bi.genome1.biased.quantile.forth = counts.allelic.mean.log2FC.degenes.bi.genome1.biased  %>% filter( KOAll_genome2_log2FC < -(log2foldchange.big) )

  ##G2 retained
  counts.allelic.mean.log2FC.degenes.bi.genome2.biased.extend = counts.allelic.mean.log2FC.degenes.bi %>%  filter( (log2FoldChange > 0.5 & get(pvaluecol) <= pvalue.cutoff.loose) )  %>% filter( AllelicComparison_log2FC >= -(mono.log2fc) ) %>% filter( KOAll_genome1_log2FC <= 0) #CK_genome1_ratio >= 0.1 & counts.allelic.mean.log2FC.degenes.bi[counts.allelic.mean.log2FC.degenes.bi$log2FoldChange > 0 & counts.allelic.mean.log2FC.degenes.bi[,pvaluecol] <= pvalue.cutoff.extend  & counts.allelic.mean.log2FC.degenes.bi$AllelicComparison_log2FC >= -(mono.log2fc),] #  & counts.allelic.mean.log2FC.degenes.bi$AllelicComparison_pvalue = pvalue.cutoff #[counts.allelic.mean.log2FC.degenes.bi.genome2.biased$log2FoldChange <= -(log2foldchange),] #allelic foldchange <= -1
  counts.allelic.mean.log2FC.degenes.bi.genome2.biased = counts.allelic.mean.log2FC.degenes.bi %>%  filter( (log2FoldChange > 0.5 & get(pvaluecol) <= pvalue.cutoff.bi) )  %>% filter( AllelicComparison_log2FC >= -(mono.log2fc) ) %>% filter(CK_genome1_ratio >= 0.1 &  KOAll_genome1_log2FC <= log2foldchange.lowest) # & AllelicComparison_pvalue >= pvalue.cutoff Status=="UP", exclude G2 expressed more than 2 fold than G1, must  #(KOAll_genome1_pvalue <= fdr & KOAll_genome1_log2FC < 0  & KOAll_genome2_log2FC * 5 > KOAll_genome1_log2FC & KOAll_genome2_log2FC > -log2foldchange.small) | 
  counts.allelic.mean.log2FC.degenes.bi.genome2.biased.logfc0.5 = counts.allelic.mean.log2FC.degenes.bi.genome2.biased  %>% filter( KOAll_genome1_log2FC < -(log2foldchange.small) )  #counts.allelic.mean.log2FC.degenes.bi.genome2.biased.logfc1[apply(counts.allelic.mean.log2FC.degenes.bi.genome2.biased.logfc1[,grep("KOAll_genome1_log2FC",colnames(counts.allelic.mean.log2FC.degenes.bi.genome2.biased.logfc1)),drop=FALSE],1,function(x){any(x<= -(log2foldchange.small))}),] #G1 is down regulated 2 fold, must
  counts.allelic.mean.log2FC.degenes.bi.genome2.biased.logfc1 = counts.allelic.mean.log2FC.degenes.bi.genome2.biased  %>% filter( KOAll_genome1_log2FC < -(log2foldchange.bi) ) #counts.allelic.mean.log2FC.degenes.bi.genome2.biased.logfc1[apply(counts.allelic.mean.log2FC.degenes.bi.genome2.biased.logfc1[,grep("KOAll_genome1_log2FC",colnames(counts.allelic.mean.log2FC.degenes.bi.genome2.biased.logfc1)),drop=FALSE],1,function(x){any(x<= -(log2foldchange))}),] #G1 is down regulated 2 fold, must
  counts.allelic.mean.log2FC.degenes.bi.genome2.biased.quantile.first = counts.allelic.mean.log2FC.degenes.bi.genome2.biased  %>% filter( KOAll_genome1_log2FC > -(log2foldchange.small) & KOAll_genome1_log2FC < 0 )
  counts.allelic.mean.log2FC.degenes.bi.genome2.biased.quantile.second = counts.allelic.mean.log2FC.degenes.bi.genome2.biased  %>% filter( KOAll_genome1_log2FC > -(log2foldchange.bi) & KOAll_genome1_log2FC <= -(log2foldchange.small) )
  counts.allelic.mean.log2FC.degenes.bi.genome2.biased.quantile.third = counts.allelic.mean.log2FC.degenes.bi.genome2.biased  %>% filter( KOAll_genome1_log2FC < -(log2foldchange.bi) )
  counts.allelic.mean.log2FC.degenes.bi.genome2.biased.quantile.thirdpart = counts.allelic.mean.log2FC.degenes.bi.genome2.biased  %>% filter( KOAll_genome1_log2FC > -(log2foldchange.big) & KOAll_genome1_log2FC <= -(log2foldchange.bi) )
  counts.allelic.mean.log2FC.degenes.bi.genome2.biased.quantile.forth = counts.allelic.mean.log2FC.degenes.bi.genome2.biased  %>% filter( KOAll_genome1_log2FC < -(log2foldchange.big) )

  counts.allelic.mean.log2FC.degenes.bi.order = rbind(counts.allelic.mean.log2FC.degenes.bi.genome1.biased, counts.allelic.mean.log2FC.degenes.bi.genome2.biased)
  
  ###########
  #Mono to none (WT_G1+1 / WT_G2+1 >= 10 (Genome1 expressed and down regulated) | WT_G2+1 / WT_G1+1 >= 10 (Genome2 expressed and down regulated))
  #mono.genome1.biased: g1_log2FC<0 & g2_log2FC < 2 (mono.genome2.biased: g2_log2FC<0 & g1_log2FC < 2)
  ###########
  counts.allelic.mean.log2FC.degenes.mono = counts.allelic.mean.log2FC.combined.degenes[!rownames(counts.allelic.mean.log2FC.combined.degenes) %in% rownames(counts.allelic.mean.log2FC.degenes.bi.order),]
  counts.allelic.mean.log2FC.degenes.mono = counts.allelic.mean.log2FC.degenes.mono[!is.na(counts.allelic.mean.log2FC.degenes.mono$AllelicComparison_pvalue),]
  
  ##G1 biased
  counts.allelic.mean.log2FC.degenes.mono.genome1.biased.extend = counts.allelic.mean.log2FC.degenes.mono %>%  filter( (log2FoldChange > 0.5 & get(pvaluecol) <= pvalue.cutoff.loose) )  %>% filter( AllelicComparison_log2FC <= -(mono.log2fc) & AllelicComparison_pvalue <= pvalue.cutoff.loose) %>% filter( KOAll_genome1_log2FC <= 0 ) #KOAll_genome1_pvalue <= fdr & 
  counts.allelic.mean.log2FC.degenes.mono.genome1.biased = counts.allelic.mean.log2FC.degenes.mono %>%  filter( (log2FoldChange > 0.5 & get(pvaluecol) <= pvalue.cutoff) )  %>% filter( AllelicComparison_log2FC <= -(mono.log2fc) & AllelicComparison_pvalue <= pvalue.cutoff) %>% filter( KOAll_genome1_log2FC <= log2foldchange.lowest ) #KOAll_genome1_pvalue <= fdr & 
  counts.allelic.mean.log2FC.degenes.mono.genome1.biased.logfc0.5 = counts.allelic.mean.log2FC.degenes.mono.genome1.biased  %>% filter( KOAll_genome1_log2FC < -(log2foldchange.small) )  #[apply(counts.allelic.mean.log2FC.degenes.mono.genome1.biased[,grep("KOAll_genome1_log2FC",colnames(counts.allelic.mean.log2FC.degenes.mono.genome1.biased)),drop=FALSE],1,function(x){all(x <= (-log2foldchange.small) )}),] #G1 is down regulated
  counts.allelic.mean.log2FC.degenes.mono.genome1.biased.logfc1 = counts.allelic.mean.log2FC.degenes.mono.genome1.biased  %>% filter( KOAll_genome1_log2FC < -(log2foldchange) )  #[apply(counts.allelic.mean.log2FC.degenes.mono.genome1.biased[,grep("KOAll_genome1_log2FC",colnames(counts.allelic.mean.log2FC.degenes.mono.genome1.biased)),drop=FALSE],1,function(x){all(x <= -log2foldchange)}),] #G1 is down regulated
  counts.allelic.mean.log2FC.degenes.mono.genome1.biased.quantile.first = counts.allelic.mean.log2FC.degenes.mono.genome1.biased  %>% filter( KOAll_genome1_log2FC > -(log2foldchange.small) & KOAll_genome1_log2FC < 0 )
  counts.allelic.mean.log2FC.degenes.mono.genome1.biased.quantile.second = counts.allelic.mean.log2FC.degenes.mono.genome1.biased  %>% filter( KOAll_genome1_log2FC > -(log2foldchange) & KOAll_genome1_log2FC <= -(log2foldchange.small) )
  counts.allelic.mean.log2FC.degenes.mono.genome1.biased.quantile.forth = counts.allelic.mean.log2FC.degenes.mono.genome1.biased  %>% filter( KOAll_genome1_log2FC < -(log2foldchange.big) & AllelicComparison_log2FC < -(log2foldchange) ) # & KOAll_genome2_log2FC > -1
  counts.allelic.mean.log2FC.degenes.mono.genome1.biased.quantile.third = counts.allelic.mean.log2FC.degenes.mono.genome1.biased  %>% filter( KOAll_genome1_log2FC < -(log2foldchange)  )
  counts.allelic.mean.log2FC.degenes.mono.genome1.biased.quantile.thirdpart = counts.allelic.mean.log2FC.degenes.mono.genome1.biased.quantile.third[!(rownames(counts.allelic.mean.log2FC.degenes.mono.genome1.biased.quantile.third) %in% rownames(counts.allelic.mean.log2FC.degenes.mono.genome1.biased.quantile.forth)), ]
  
  ##G2 biased
  counts.allelic.mean.log2FC.degenes.mono.genome2.biased.extend = counts.allelic.mean.log2FC.degenes.mono %>%  filter( (log2FoldChange < -0.5 & get(pvaluecol) <= pvalue.cutoff.loose) )  %>% filter( AllelicComparison_log2FC >= (mono.log2fc) & AllelicComparison_pvalue <= pvalue.cutoff.loose) %>% filter( KOAll_genome2_log2FC <= 0) #KOAll_genome2_pvalue <= fdr & 
  counts.allelic.mean.log2FC.degenes.mono.genome2.biased = counts.allelic.mean.log2FC.degenes.mono %>%  filter( (log2FoldChange < -0.5 & get(pvaluecol) <= pvalue.cutoff) )  %>% filter( AllelicComparison_log2FC >= (mono.log2fc) & AllelicComparison_pvalue <= pvalue.cutoff) %>% filter( KOAll_genome2_log2FC <= log2foldchange.lowest ) #KOAll_genome2_pvalue <= pvalue.cutoff & 
  counts.allelic.mean.log2FC.degenes.mono.genome2.biased.logfc0.5 = counts.allelic.mean.log2FC.degenes.mono.genome2.biased %>% filter( KOAll_genome2_log2FC < -(log2foldchange.small) )  #[apply(counts.allelic.mean.log2FC.degenes.mono.genome2.biased[,grep("KOAll_genome2_log2FC",colnames(counts.allelic.mean.log2FC.degenes.mono.genome2.biased)),drop=FALSE],1,function(x){all(x <= (-log2foldchange.small) )}),] #G2 is down regulated
  counts.allelic.mean.log2FC.degenes.mono.genome2.biased.logfc1 = counts.allelic.mean.log2FC.degenes.mono.genome2.biased %>% filter( KOAll_genome2_log2FC < -(log2foldchange) )  #[apply(counts.allelic.mean.log2FC.degenes.mono.genome2.biased[,grep("KOAll_genome2_log2FC",colnames(counts.allelic.mean.log2FC.degenes.mono.genome2.biased)),drop=FALSE],1,function(x){all(x <= -log2foldchange)}),] #G2 is down regulated
  counts.allelic.mean.log2FC.degenes.mono.genome2.biased.quantile.first = counts.allelic.mean.log2FC.degenes.mono.genome2.biased  %>% filter( KOAll_genome2_log2FC > -(log2foldchange.small) & KOAll_genome2_log2FC < 0 )
  counts.allelic.mean.log2FC.degenes.mono.genome2.biased.quantile.second = counts.allelic.mean.log2FC.degenes.mono.genome2.biased  %>% filter( KOAll_genome2_log2FC > -(log2foldchange) & KOAll_genome2_log2FC <= -(log2foldchange.small) )
  counts.allelic.mean.log2FC.degenes.mono.genome2.biased.quantile.forth = counts.allelic.mean.log2FC.degenes.mono.genome2.biased  %>% filter( KOAll_genome2_log2FC < -(log2foldchange.big)  & AllelicComparison_log2FC > (log2foldchange) ) # & KOAll_genome1_log2FC > -1
  counts.allelic.mean.log2FC.degenes.mono.genome2.biased.quantile.third = counts.allelic.mean.log2FC.degenes.mono.genome2.biased  %>% filter( KOAll_genome2_log2FC < -(log2foldchange)) # & CK_genome1 > expression.cutoff 
  counts.allelic.mean.log2FC.degenes.mono.genome2.biased.quantile.thirdpart = counts.allelic.mean.log2FC.degenes.mono.genome2.biased.quantile.third[!(rownames(counts.allelic.mean.log2FC.degenes.mono.genome2.biased.quantile.third) %in% rownames(counts.allelic.mean.log2FC.degenes.mono.genome2.biased.quantile.forth)), ]
  
  counts.allelic.mean.log2FC.degenes.mono.order = rbind(counts.allelic.mean.log2FC.degenes.mono.genome2.biased, counts.allelic.mean.log2FC.degenes.mono.genome1.biased)
  
  
  ###########
  ##Bi to none (changed the first two lines, not run yet)
  ###########
  counts.allelic.mean.log2FC.degenes.bi.both.down = counts.allelic.mean.log2FC.combined.degenes[!rownames(counts.allelic.mean.log2FC.combined.degenes) %in% c(rownames(counts.allelic.mean.log2FC.degenes.bi.order),rownames(counts.allelic.mean.log2FC.degenes.mono.order)),]
  counts.allelic.mean.log2FC.degenes.bi.both.down.extend = counts.allelic.mean.log2FC.degenes.bi.both.down %>%  filter(basic_Status=="DOWN" & get(pvaluecol) > pvalue.cutoff.loose) %>% filter(KOAll_genome1_pvalue <= pvalue.cutoff.loose | KOAll_genome2_pvalue <= pvalue.cutoff.loose) %>% filter(KOAll_genome1_log2FC <= 0 | KOAll_genome2_log2FC <= 0)
  counts.allelic.mean.log2FC.degenes.bi.both.down = counts.allelic.mean.log2FC.degenes.bi.both.down %>%  filter(basic_Status=="DOWN" & get(pvaluecol) > pvalue.cutoff) %>% filter(KOAll_genome1_pvalue <= pvalue.cutoff | KOAll_genome2_pvalue <= pvalue.cutoff) %>% filter(KOAll_genome1_log2FC <= log2foldchange.lowest | KOAll_genome2_log2FC <= log2foldchange.lowest)
  counts.allelic.mean.log2FC.degenes.bi.both.down.quantile.first = counts.allelic.mean.log2FC.degenes.bi.both.down %>% filter(KOAll_genome1_log2FC > -(log2foldchange.small) | KOAll_genome2_log2FC > -(log2foldchange.small) )
  counts.allelic.mean.log2FC.degenes.bi.both.down.quantile.second = counts.allelic.mean.log2FC.degenes.bi.both.down %>% filter((KOAll_genome1_log2FC <= -(log2foldchange.small) & KOAll_genome1_log2FC > -(log2foldchange)) | (KOAll_genome2_log2FC <= -(log2foldchange.small) & KOAll_genome2_log2FC > -(log2foldchange)) )
  counts.allelic.mean.log2FC.degenes.bi.both.down.quantile.second = counts.allelic.mean.log2FC.degenes.bi.both.down.quantile.second[!(rownames(counts.allelic.mean.log2FC.degenes.bi.both.down.quantile.second) %in% rownames(counts.allelic.mean.log2FC.degenes.bi.both.down.quantile.first)), ] 
  counts.allelic.mean.log2FC.degenes.bi.both.down.quantile.forth = counts.allelic.mean.log2FC.degenes.bi.both.down %>% filter(KOAll_genome1_log2FC <= -(log2foldchange.big) & KOAll_genome2_log2FC <= -(log2foldchange.big) )
  counts.allelic.mean.log2FC.degenes.bi.both.down.quantile.third = counts.allelic.mean.log2FC.degenes.bi.both.down[! rownames(counts.allelic.mean.log2FC.degenes.bi.both.down) %in% rownames(rbind(counts.allelic.mean.log2FC.degenes.bi.both.down.quantile.first, counts.allelic.mean.log2FC.degenes.bi.both.down.quantile.second)),]
  counts.allelic.mean.log2FC.degenes.bi.both.down.quantile.thirdpart = counts.allelic.mean.log2FC.degenes.bi.both.down[! rownames(counts.allelic.mean.log2FC.degenes.bi.both.down) %in% rownames(rbind(counts.allelic.mean.log2FC.degenes.bi.both.down.quantile.first, counts.allelic.mean.log2FC.degenes.bi.both.down.quantile.second, counts.allelic.mean.log2FC.degenes.bi.both.down.quantile.forth)),]
  counts.allelic.mean.log2FC.degenes.bi.both.down.logfc0.5 = as.data.frame(rbind(counts.allelic.mean.log2FC.degenes.bi.both.down.quantile.second, counts.allelic.mean.log2FC.degenes.bi.both.down.quantile.third))
  counts.allelic.mean.log2FC.degenes.bi.both.down.logfc1 = counts.allelic.mean.log2FC.degenes.bi.both.down.quantile.third

  ##export to bed
  export(species.database.anno.gene.granges[rownames(counts.allelic.mean.log2FC.degenes.mono.genome1.biased)],con=file.path(output.dir,"DE.allelic.mono.genome1.biased.all.bed"))
  export(species.database.anno.gene.granges[rownames(counts.allelic.mean.log2FC.degenes.mono.genome1.biased.logfc0.5)],con=file.path(output.dir,"DE.allelic.mono.genome1.biased.bed"))
  export(species.database.anno.gene.granges[rownames(counts.allelic.mean.log2FC.degenes.mono.genome1.biased.logfc1)],con=file.path(output.dir,"DE.allelic.mono.genome1.biased.logfc1.bed"))
  export(species.database.anno.gene.granges[rownames(counts.allelic.mean.log2FC.degenes.mono.genome2.biased)],con=file.path(output.dir,"DE.allelic.mono.genome2.biased.all.bed"))
  export(species.database.anno.gene.granges[rownames(counts.allelic.mean.log2FC.degenes.mono.genome2.biased.logfc0.5)],con=file.path(output.dir,"DE.allelic.mono.genome2.biased.bed"))
  export(species.database.anno.gene.granges[rownames(counts.allelic.mean.log2FC.degenes.mono.genome2.biased.logfc1)],con=file.path(output.dir,"DE.allelic.mono.genome2.biased.logfc1.bed"))
  export(species.database.anno.gene.granges[rownames(counts.allelic.mean.log2FC.degenes.bi.genome1.biased)],con=file.path(output.dir,"DE.allelic.bi.genome1.biased.all.bed"))
  export(species.database.anno.gene.granges[rownames(counts.allelic.mean.log2FC.degenes.bi.genome1.biased.logfc0.5)],con=file.path(output.dir,"DE.allelic.bi.genome1.biased.bed"))
  export(species.database.anno.gene.granges[rownames(counts.allelic.mean.log2FC.degenes.bi.genome1.biased.logfc1)],con=file.path(output.dir,"DE.allelic.bi.genome1.biased.logfc1.bed"))
  export(species.database.anno.gene.granges[rownames(counts.allelic.mean.log2FC.degenes.bi.genome2.biased)],con=file.path(output.dir,"DE.allelic.bi.genome2.biased.all.bed"))
  export(species.database.anno.gene.granges[rownames(counts.allelic.mean.log2FC.degenes.bi.genome2.biased.logfc0.5)],con=file.path(output.dir,"DE.allelic.bi.genome2.biased.bed"))
  export(species.database.anno.gene.granges[rownames(counts.allelic.mean.log2FC.degenes.bi.genome2.biased.logfc1)],con=file.path(output.dir,"DE.allelic.bi.genome2.biased.logfc1.bed"))
  export(species.database.anno.gene.granges[rownames(counts.allelic.mean.log2FC.degenes.bi.both.down)],con=file.path(output.dir,"DE.allelic.bi.both.down.all.bed"))
  export(species.database.anno.gene.granges[rownames(counts.allelic.mean.log2FC.degenes.bi.both.down.logfc0.5)],con=file.path(output.dir,"DE.allelic.bi.both.down.bed"))
  export(species.database.anno.gene.granges[rownames(counts.allelic.mean.log2FC.degenes.bi.both.down.logfc1)],con=file.path(output.dir,"DE.allelic.bi.both.down.logfc1.bed"))
  
  build.export.quantile = function(category="bi.genome2.biased"){
    for (filename in c("first", "second", "thirdpart", "forth")){
      export(species.database.anno.gene.granges[rownames(get(paste0("counts.allelic.mean.log2FC.degenes.", category, ".quantile.", filename)))],con=file.path(output.dir,paste0("quantile.", gsub("part", "", filename), ".DE.allelic.", category, ".bed")))
    }
  }
  build.export.quantile(category="bi.genome2.biased")
  build.export.quantile(category="bi.genome1.biased")
  build.export.quantile(category="mono.genome2.biased")
  build.export.quantile(category="mono.genome1.biased")
  build.export.quantile(category="bi.both.down")
  
  
  ##combine all the categories and store the table
  counts.allelic.mean.log2FC.combined.degenes.all = counts.allelic.mean.log2FC.combined.degenes[c(rownames(counts.allelic.mean.log2FC.degenes.bi.genome2.biased),
                                                                                              rownames(counts.allelic.mean.log2FC.degenes.bi.genome1.biased),
                                                                                              rownames(counts.allelic.mean.log2FC.degenes.mono.genome1.biased),
                                                                                              rownames(counts.allelic.mean.log2FC.degenes.mono.genome2.biased),
                                                                                              rownames(counts.allelic.mean.log2FC.degenes.bi.both.down)), ]
  counts.allelic.mean.log2FC.combined.degenes.all$Quantile = "log2FC > -1"
  counts.allelic.mean.log2FC.combined.degenes.all$AllelicStatus = c(rep("Bi_Allele2_biased",nrow(counts.allelic.mean.log2FC.degenes.bi.genome2.biased)),
                                                                rep("Bi_Allele1_biased",nrow(counts.allelic.mean.log2FC.degenes.bi.genome1.biased)),
                                                                rep("Mono_Allele1_biased",nrow(counts.allelic.mean.log2FC.degenes.mono.genome1.biased)),
                                                                rep("Mono_Allele2_biased",nrow(counts.allelic.mean.log2FC.degenes.mono.genome2.biased)),
                                                                rep("Bi_both_down",nrow(counts.allelic.mean.log2FC.degenes.bi.both.down)))
  write.table(counts.allelic.mean.log2FC.combined.degenes.all,file = file.path(output.dir,"Expression","log2FC.separated.genome.degenes.all.txt"), sep="\t",quote = FALSE, col.names = NA)
  print("DEgene: All DEgenes in the log2FC.separated.genome.degenes table, pvalue < 0.01, and mono.logfc and genome1 percentage")
  print(nrow(counts.allelic.mean.log2FC.combined.degenes.all))
  log2FC.separated.genome.degenes.all = counts.allelic.mean.log2FC.combined.degenes.all
  save(log2FC.separated.genome.degenes.all,file = file.path(output.dir,"Expression","log2FC.separated.genome.degenes.all.RData"))
  print(table(log2FC.separated.genome.degenes.all$AllelicStatus))
  
  ##combine all the categories that the logfc>0.5 and store the table
  counts.allelic.mean.log2FC.combined.degenes.logfc0.5 = counts.allelic.mean.log2FC.combined.degenes[c(rownames(counts.allelic.mean.log2FC.degenes.bi.genome2.biased.logfc0.5),
                                                                                                       rownames(counts.allelic.mean.log2FC.degenes.bi.genome1.biased.logfc0.5),
                                                                                                       rownames(counts.allelic.mean.log2FC.degenes.mono.genome1.biased.logfc0.5),
                                                                                                       rownames(counts.allelic.mean.log2FC.degenes.mono.genome2.biased.logfc0.5),
                                                                                                       rownames(counts.allelic.mean.log2FC.degenes.bi.both.down.logfc0.5)), ]
  counts.allelic.mean.log2FC.combined.degenes.logfc0.5$AllelicStatus = c(rep("Bi_Allele2_biased",nrow(counts.allelic.mean.log2FC.degenes.bi.genome2.biased.logfc0.5)),
                                                                         rep("Bi_Allele1_biased",nrow(counts.allelic.mean.log2FC.degenes.bi.genome1.biased.logfc0.5)),
                                                                         rep("Mono_Allele1_biased",nrow(counts.allelic.mean.log2FC.degenes.mono.genome1.biased.logfc0.5)),
                                                                         rep("Mono_Allele2_biased",nrow(counts.allelic.mean.log2FC.degenes.mono.genome2.biased.logfc0.5)),
                                                                         rep("Bi_both_down",nrow(counts.allelic.mean.log2FC.degenes.bi.both.down.logfc0.5)) )
  write.table(counts.allelic.mean.log2FC.combined.degenes.logfc0.5,file = file.path(output.dir,"Expression","log2FC.separated.genome.degenes.txt"), sep="\t",quote = FALSE, col.names = NA)
  print("DEgene: log2FC > 0.5 , pvalue < 0.01 DEgenes in the log2FC.separated.genome.degenes.logfc0.5 table")
  print(nrow(counts.allelic.mean.log2FC.combined.degenes.logfc0.5))
  log2FC.separated.genome.degenes = counts.allelic.mean.log2FC.combined.degenes.logfc0.5
  save(log2FC.separated.genome.degenes,file = file.path(output.dir,"Expression","log2FC.separated.genome.degenes.RData"))
  print(table(log2FC.separated.genome.degenes$AllelicStatus))
  
  ##combine all the categories that the logfc>1 and store the table
  counts.allelic.mean.log2FC.combined.degenes.logfc1 = counts.allelic.mean.log2FC.combined.degenes[c(rownames(counts.allelic.mean.log2FC.degenes.bi.genome2.biased.logfc1),
                                                                                                     rownames(counts.allelic.mean.log2FC.degenes.bi.genome1.biased.logfc1),
                                                                                                     rownames(counts.allelic.mean.log2FC.degenes.mono.genome1.biased.logfc1),
                                                                                                     rownames(counts.allelic.mean.log2FC.degenes.mono.genome2.biased.logfc1),
                                                                                                     rownames(counts.allelic.mean.log2FC.degenes.bi.both.down.logfc1)), ]
  counts.allelic.mean.log2FC.combined.degenes.logfc1$Quantile = "log2FC > -1"
  counts.allelic.mean.log2FC.combined.degenes.logfc1$AllelicStatus = c(rep("Bi_Allele2_biased",nrow(counts.allelic.mean.log2FC.degenes.bi.genome2.biased.logfc1)),
                                                                       rep("Bi_Allele1_biased",nrow(counts.allelic.mean.log2FC.degenes.bi.genome1.biased.logfc1)),
                                                                       rep("Mono_Allele1_biased",nrow(counts.allelic.mean.log2FC.degenes.mono.genome1.biased.logfc1)),
                                                                       rep("Mono_Allele2_biased",nrow(counts.allelic.mean.log2FC.degenes.mono.genome2.biased.logfc1)),
                                                                       rep("Bi_both_down",nrow(counts.allelic.mean.log2FC.degenes.bi.both.down.logfc1)) )
  write.table(counts.allelic.mean.log2FC.combined.degenes.logfc1,file = file.path(output.dir,"Expression","log2FC.separated.genome.degenes.logfc1.txt"), sep="\t",quote = FALSE, col.names = NA)
  print("DEgene: log2FC > 1, FDR < 0.1 DEgenes in the log2FC.separated.genome.degenes.logfc1 table")
  print(nrow(counts.allelic.mean.log2FC.combined.degenes.logfc1))
  # export(species.database.anno.gene.granges[rownames(counts.allelic.mean.log2FC.combined.degenes)],con=file.path(output.dir,"DE.allelic.all.degenes.bed"))
  log2FC.separated.genome.degenes.logfc1 = counts.allelic.mean.log2FC.combined.degenes.logfc1
  save(log2FC.separated.genome.degenes.logfc1,file = file.path(output.dir,"Expression","log2FC.separated.genome.degenes.logfc1.RData"))
  print(table(log2FC.separated.genome.degenes.logfc1$AllelicStatus))
  
  ##combine all the categories that the P < 0.1 and store the table
  counts.allelic.mean.log2FC.combined.degenes.extend = counts.allelic.mean.log2FC.combined.degenes[c(rownames(counts.allelic.mean.log2FC.degenes.bi.genome2.biased.extend),
                                                                                                   rownames(counts.allelic.mean.log2FC.degenes.bi.genome1.biased.extend),
                                                                                                   rownames(counts.allelic.mean.log2FC.degenes.mono.genome1.biased.extend),
                                                                                                   rownames(counts.allelic.mean.log2FC.degenes.mono.genome2.biased.extend),
                                                                                                   rownames(counts.allelic.mean.log2FC.degenes.bi.both.down.extend)), ]
  counts.allelic.mean.log2FC.combined.degenes.extend$Quantile = "log2FC > -1"
  counts.allelic.mean.log2FC.combined.degenes.extend$AllelicStatus = c(rep("Bi_Allele2_biased",nrow(counts.allelic.mean.log2FC.degenes.bi.genome2.biased.extend)),
                                                                     rep("Bi_Allele1_biased",nrow(counts.allelic.mean.log2FC.degenes.bi.genome1.biased.extend)),
                                                                     rep("Mono_Allele1_biased",nrow(counts.allelic.mean.log2FC.degenes.mono.genome1.biased.extend)),
                                                                     rep("Mono_Allele2_biased",nrow(counts.allelic.mean.log2FC.degenes.mono.genome2.biased.extend)),
                                                                     rep("Bi_both_down",nrow(counts.allelic.mean.log2FC.degenes.bi.both.down.extend)) )
  write.table(counts.allelic.mean.log2FC.combined.degenes.extend,file = file.path(output.dir,"Expression","log2FC.separated.genome.degenes.extend.txt"), sep="\t",quote = FALSE, col.names = NA)
  print("DEgene: All DEgenes in the log2FC.separated.genome.degenes table, pvalue < 0.1, and mono.logfc and genome1 percentage")
  print(nrow(counts.allelic.mean.log2FC.combined.degenes.extend))
  # export(species.database.anno.gene.granges[rownames(counts.allelic.mean.log2FC.combined.degenes)],con=file.path(output.dir,"DE.allelic.all.degenes.bed"))
  log2FC.separated.genome.degenes.extend = counts.allelic.mean.log2FC.combined.degenes.extend
  save(log2FC.separated.genome.degenes.extend,file = file.path(output.dir,"Expression","log2FC.separated.genome.degenes.extend.RData"))
  print(table(log2FC.separated.genome.degenes.extend$AllelicStatus))
  
  ##combine all the categories and store the table
  merged.rownames = c()
  category.column.name = c()
  filename.column.name = c()
  for (category in c("bi.genome1.biased", "bi.genome2.biased", "mono.genome1.biased", "mono.genome2.biased", "bi.both.down")){
    for (filename in c("first", "second", "third")){ #, "forth"
      merged.rownames = c(merged.rownames, rownames(get(paste0("counts.allelic.mean.log2FC.degenes.", category, ".quantile.", filename))) )
      category.column.name = c(category.column.name, rep(category,nrow(get(paste0("counts.allelic.mean.log2FC.degenes.", category, ".quantile.", filename)))))
      filename.column.name = c(filename.column.name, rep(toupper(filename),nrow(get(paste0("counts.allelic.mean.log2FC.degenes.", category, ".quantile.", filename)))))
      filename.column.name = gsub("FIRST", "log2FC > -1", gsub("SECOND", "-2< log2FC < -1", gsub("THIRD", "log2FC < -2", filename.column.name)))  #gsub("THIRD", "2< log2FC < 3", gsub("FORTH", "log2FC > 3", 
    } 
  }
  counts.allelic.mean.log2FC.combined.degenes.quantile = counts.allelic.mean.log2FC.combined.degenes[merged.rownames, ]
  counts.allelic.mean.log2FC.combined.degenes.quantile$Quantile = filename.column.name
  counts.allelic.mean.log2FC.combined.degenes.quantile$AllelicStatus = gsub("[.]", "_", gsub("genome", "Allele", gsub("^bi", "Bi", gsub("mono", "Mono", gsub("bi.both.down", "Bi_both_down", category.column.name)))))
  write.table(counts.allelic.mean.log2FC.combined.degenes.quantile,file = file.path(output.dir,"Expression","log2FC.separated.genome.degenes.quantile.txt"), sep="\t",quote = FALSE, col.names = NA)
  print("DEgene: Quantile 3 DEgenes in the log2FC.separated.genome.degenes table, pvalue < 0.05, and mono.logfc and genome1 percentage")
  print(nrow(counts.allelic.mean.log2FC.combined.degenes.quantile))
  # export(species.database.anno.gene.granges[rownames(counts.allelic.mean.log2FC.combined.degenes)],con=file.path(output.dir,"DE.allelic.quantile.degenes.bed"))
  log2FC.separated.genome.degenes.quantile = counts.allelic.mean.log2FC.combined.degenes.quantile
  save(log2FC.separated.genome.degenes.quantile,file = file.path(output.dir,"Expression","log2FC.separated.genome.degenes.quantile.RData"))
  print(table(log2FC.separated.genome.degenes.quantile$AllelicStatus, log2FC.separated.genome.degenes.quantile$Quantile))
  
  ##combine quantilepart the categories and store the table
  merged.rownames = c()
  category.column.name = c()
  filename.column.name = c()
  for (category in c("bi.genome1.biased", "bi.genome2.biased", "mono.genome1.biased", "mono.genome2.biased", "bi.both.down")){
    for (filename in c("first", "second", "thirdpart", "forth")){ #
      merged.rownames = c(merged.rownames, rownames(get(paste0("counts.allelic.mean.log2FC.degenes.", category, ".quantile.", filename))) )
      category.column.name = c(category.column.name, rep(category,nrow(get(paste0("counts.allelic.mean.log2FC.degenes.", category, ".quantile.", filename)))))
      filename.column.name = c(filename.column.name, rep(toupper(filename),nrow(get(paste0("counts.allelic.mean.log2FC.degenes.", category, ".quantile.", filename)))))
      filename.column.name = gsub("FIRST", "log2FC > -1", gsub("SECOND", "-2< log2FC < -1", gsub("THIRDPART", "-3< log2FC < -2", gsub("FORTH", "log2FC < -3",filename.column.name))))  # 
    } 
  }
  counts.allelic.mean.log2FC.combined.degenes.quantile = counts.allelic.mean.log2FC.combined.degenes[merged.rownames, ]
  counts.allelic.mean.log2FC.combined.degenes.quantile$Quantile = filename.column.name
  counts.allelic.mean.log2FC.combined.degenes.quantile$AllelicStatus = gsub("[.]", "_", gsub("genome", "Allele", gsub("^bi", "Bi", gsub("mono", "Mono", gsub("bi.both.down", "Bi_both_down", category.column.name)))))
  write.table(counts.allelic.mean.log2FC.combined.degenes.quantile,file = file.path(output.dir,"Expression","log2FC.separated.genome.degenes.quantilepart.txt"), sep="\t",quote = FALSE, col.names = NA)
  print("DEgene: Quantile 4 DEgenes in the log2FC.separated.genome.degenes table, pvalue < 0.05, and mono.logfc and genome1 percentage")
  print(nrow(counts.allelic.mean.log2FC.combined.degenes.quantile))
  # export(species.database.anno.gene.granges[rownames(counts.allelic.mean.log2FC.combined.degenes)],con=file.path(output.dir,"DE.allelic.all.degenes.bed"))
  log2FC.separated.genome.degenes.quantilepart = counts.allelic.mean.log2FC.combined.degenes.quantile
  save(log2FC.separated.genome.degenes.quantilepart,file = file.path(output.dir,"Expression","log2FC.separated.genome.degenes.quantilepart.RData"))
  print(table(log2FC.separated.genome.degenes.quantilepart$AllelicStatus, log2FC.separated.genome.degenes.quantilepart$Quantile))
  
  ##All mono allelic genes
  mono.allelic.genes = read.delim(file.path(prewd, deseq2.directory,"DEseq_resutls_CK_two_genomes.tsv"),row.names = 1)
  gene.mono.g1.biased = mono.allelic.genes[mono.allelic.genes$Status=="DOWN",]
  gene.mono.g2.biased = mono.allelic.genes[mono#"/data/akhtar/Mouse2019AlleleSpecific/NPCBlCa2019RNAseqMouse".allelic.genes$Status=="UP",]
  
  export(species.database.anno.gene.granges[names(species.database.anno.gene.granges) %in% rownames(gene.mono.g1.biased)],con=file.path(output.dir,"DE.all.mono.genome1.biased.bed"))
  export(species.database.anno.gene.granges[names(species.database.anno.gene.granges) %in% rownames(gene.mono.g2.biased)],con=file.path(output.dir,"DE.all.mono.genome2.biased.bed"))
  
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}