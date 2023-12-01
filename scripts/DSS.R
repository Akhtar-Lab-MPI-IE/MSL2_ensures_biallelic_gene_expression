require("bsseq")
require("BiocParallel")
require("DSS")
require("ggplot2")

#Set threads
BPPARAM = MulticoreParam(workers = 16, progressbar=TRUE)

#Set workdir
setwd('methExtract')

#Read up the sample sheet.
ss <- read.delim("sampleSheet.tsv")

#Define cytosine reports.
CpGreports <- c("KO1_6_NPCs_Rep1_R1_bismark_bt2_pe.Cast.CpG_report.txt",
                "KO1_6_NPCs_Rep1_R1_bismark_bt2_pe.Svj.CpG_report.txt",
                "KO1_6_NPCs_Rep2_R1_bismark_bt2_pe.Cast.CpG_report.txt",
                "KO1_6_NPCs_Rep2_R1_bismark_bt2_pe.Svj.CpG_report.txt",
                "KO2_16_NPCs_Rep1_R1_bismark_bt2_pe.Cast.CpG_report.txt",
                "KO2_16_NPCs_Rep1_R1_bismark_bt2_pe.Svj.CpG_report.txt",
                "KO2_16_NPCs_Rep2_R1_bismark_bt2_pe.Cast.CpG_report.txt",
                "KO2_16_NPCs_Rep2_R1_bismark_bt2_pe.Svj.CpG_report.txt",
                "WT_G4_NPCs_Rep1_R1_bismark_bt2_pe.Cast.CpG_report.txt",
                "WT_G4_NPCs_Rep1_R1_bismark_bt2_pe.Svj.CpG_report.txt",
                "WT_G4_NPCs_Rep2_R1_bismark_bt2_pe.Cast.CpG_report.txt",
                "WT_G4_NPCs_Rep2_R1_bismark_bt2_pe.Svj.CpG_report.txt",
                "WT_G4_NPCs_Rep3_R1_bismark_bt2_pe.Cast.CpG_report.txt",
                "WT_G4_NPCs_Rep3_R1_bismark_bt2_pe.Svj.CpG_report.txt",
                "WT_G4_NPCs_Rep4_R1_bismark_bt2_pe.Cast.CpG_report.txt",
                "WT_G4_NPCs_Rep4_R1_bismark_bt2_pe.Svj.CpG_report.txt")

#Read the files.
bs <- read.bismark(CpGreports,
                   BPPARAM=BPPARAM,
                  colData = ss)

#Incorporate samplenames
sampleNames(bs) <- gsub(".CpG_report.txt","", gsub("_R1_bismark_bt2_pe", "", as.character(ss$name)))

# Knockouts - Svj (group1) vs Cast (group2)

kobs_dmlT = DMLtest(bs,
                    group1=c("KO1_6_NPCs_Rep1.Svj",
                             "KO1_6_NPCs_Rep2.Svj",
                             "KO2_16_NPCs_Rep1.Svj",
                             "KO2_16_NPCs_Rep2.Svj"),
                    group2=c("KO1_6_NPCs_Rep1.Cast",
                             "KO1_6_NPCs_Rep2.Cast",
                             "KO2_16_NPCs_Rep1.Cast",
                             "KO2_16_NPCs_Rep2.Cast"),
                    smoothing=TRUE,
                    BPPARAM = BPPARAM)

kobs_dmls <- callDML(kobs_dmlT)
kobs_dmrs <- callDMR(kobs_dmlT, delta=0.1, dis.merge=300, minCG=10)


wtbs_dmlT = DMLtest(bs,
                    group1=c("WT_G4_NPCs_Rep1.Svj",
                             "WT_G4_NPCs_Rep2.Svj",
                             "WT_G4_NPCs_Rep3.Svj",
                             "WT_G4_NPCs_Rep4.Svj"),
                    group2=c("WT_G4_NPCs_Rep1.Cast",
                             "WT_G4_NPCs_Rep2.Cast",
                             "WT_G4_NPCs_Rep3.Cast",
                             "WT_G4_NPCs_Rep4.Cast"),
                    smoothing=TRUE,
                    BPPARAM = BPPARAM)

wtbs_dmls <- callDML(wtbs_dmlT)
wtbs_dmrs <- callDMR(wtbs_dmlT, delta=0.1, dis.merge=300, minCG=10)

# Write out data
## DMLs
write.csv(kobs_dmls, 'ko_svjvscast_dmls.csv', quote=F)
write.csv(wtbs_dmls, 'wt_svjvscast_dmls.csv', quote=F)
## DMRs
write.csv(kobs_dmrs, 'ko_svjvscast_dmrs.csv', quote=F)
write.csv(wtbs_dmrs, 'wt_svjvscast_dmrs.csv', quote=F)





CastKOvsWT = DMLtest(bs,
                    group1=c("KO1_6_NPCs_Rep1.Cast",
                             "KO1_6_NPCs_Rep2.Cast",
                             "KO2_16_NPCs_Rep1.Cast",
                             "KO2_16_NPCs_Rep2.Cast"),
                    group2=c("WT_G4_NPCs_Rep1.Cast",
                             "WT_G4_NPCs_Rep2.Cast",
                             "WT_G4_NPCs_Rep3.Cast",
                             "WT_G4_NPCs_Rep4.Cast"),
                    smoothing=TRUE,
                    BPPARAM = BPPARAM)

CastKOvsWT_dmls <- callDML(CastKOvsWT)
CastKOvsWT_dmrs <- callDMR(CastKOvsWT, delta=0.1, dis.merge=300, minCG=10)
write.csv(CastKOvsWT_dmls, 'CastKOvsWT_dmls.csv', quote=F)
write.csv(CastKOvsWT_dmrs, 'CastKOvsWT_dmrs.csv', quote=F)

SvjKOvsWT = DMLtest(bs,
                     group1=c("KO1_6_NPCs_Rep1.Svj",
                              "KO1_6_NPCs_Rep2.Svj",
                              "KO2_16_NPCs_Rep1.Svj",
                              "KO2_16_NPCs_Rep2.Svj"),
                     group2=c("WT_G4_NPCs_Rep1.Svj",
                              "WT_G4_NPCs_Rep2.Svj",
                              "WT_G4_NPCs_Rep3.Svj",
                              "WT_G4_NPCs_Rep4.Svj"),
                     smoothing=TRUE,
                    BPPARAM = BPPARAM)

SvjKOvsWT_dmls <- callDML(SvjKOvsWT)
SvjKOvsWT_dmrs <- callDMR(SvjKOvsWT, delta=0.1, dis.merge=300, minCG=10)
write.csv(SvjKOvsWT_dmls, 'SvjKOvsWT_dmls.csv', quote=F)
write.csv(SvjKOvsWT_dmrs, 'SvjKOvsWT_dmrs.csv', quote=F)



#Set threads
BPPARAM = MulticoreParam(workers = 16, progressbar=TRUE)

#Set workdir
setwd('../fullmethExtract')

#Read up the sample sheet.
ss <- read.delim("sampleSheet_full.tsv")

#Define cytosine reports.
CpGreports <- c("WT_G4_NPCs_Rep1_R1_bismark_bt2_pe.CpG_report.txt",
                "WT_G4_NPCs_Rep2_R1_bismark_bt2_pe.CpG_report.txt",
                "WT_G4_NPCs_Rep3_R1_bismark_bt2_pe.CpG_report.txt",
                "WT_G4_NPCs_Rep4_R1_bismark_bt2_pe.CpG_report.txt",
                "KO1_6_NPCs_Rep1_R1_bismark_bt2_pe.CpG_report.txt",
                "KO1_6_NPCs_Rep2_R1_bismark_bt2_pe.CpG_report.txt",
                "KO2_16_NPCs_Rep1_R1_bismark_bt2_pe.CpG_report.txt",
                "KO2_16_NPCs_Rep2_R1_bismark_bt2_pe.CpG_report.txt")

#Read the files.
bs <- read.bismark(CpGreports,
                   BPPARAM=BPPARAM,
                   colData = ss)

#Incorporate samplenames
sampleNames(bs) <- gsub(".CpG_report.txt","", gsub("_R1_bismark_bt2_pe", "", as.character(ss$name)))

# Knockouts - KO (group1) vs WT (group2).

KOvsWT_dmlT = DMLtest(bs,
                    group1=c("KO1_6_NPCs_Rep1",
                             "KO1_6_NPCs_Rep2",
                             "KO2_16_NPCs_Rep1",
                             "KO2_16_NPCs_Rep2"),                      
                    group2=c("WT_G4_NPCs_Rep1",
                             "WT_G4_NPCs_Rep2",
                             "WT_G4_NPCs_Rep3",
                             "WT_G4_NPCs_Rep4"),
                    smoothing=TRUE,
                    BPPARAM = BPPARAM)

KOvsWT_dmls <- callDML(KOvsWT_dmlT)
KOvsWT_dmrs <- callDMR(KOvsWT_dmlT, delta=0.1, dis.merge=300, minCG=10)


# Write out data
## DMLs
write.csv(KOvsWT_dmls, 'KOvsWT_full_dmls.csv', quote=F)
## DMRs
write.csv(KOvsWT_dmrs, 'KOvsWT_full_dmrs.csv', quote=F)
