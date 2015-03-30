rm(list=ls())

##
## Load libraries
##
require(doMC)
registerDoMC(2)

require(sqldf)
require(gdata)
require(gtools)
require(annotate)
require(snpStats)
require(data.table)
require(MatrixEQTL)
require(org.Hs.eg.db)
require(GenomicRanges)
require(VariantAnnotation)

##
## Define custom functions
##
my.t.test.p.value <- function(...) { 
  obj<-try(t.test(...), silent=TRUE) 
  if (is(obj, "try-error")) return(NA) else if (is.nan(obj$p.value)) return(NA) else return(obj$p.value) 
} 

my.wilcox.test.p.value <- function(...) { 
  obj<-try(wilcox.test(...), silent=TRUE) 
  if (is(obj, "try-error")) return(NA) else if (is.nan(obj$p.value)) return(NA) else return(obj$p.value) 
} 

##
## Initialize global variables
##
if (!is.na(pmatch("darwin", version$os))) {
    rootDir = "/Volumes/orchestra"
} else {
    rootDir = "/home/sl279"
}

chrs = c(1:22, "X", "Y")
chrs = factor(chrs, level=chrs)
cchrs = sapply(chrs, function(x) paste("chr", x, sep=""))
cchrs = factor(cchrs, level=cchrs)

baseDir = file.path(rootDir, "BiO/Research/Recuronco")
vcfDir = file.path(baseDir, "vcf")
panVcfDir = file.path(vcfDir, "PANCAN")
gdacDir = file.path(baseDir, "gdac")

gdacStdDataDate = "2014_10_17"
gdacStdDataDateNoUb = gsub("_", "", gdacStdDataDate)
gdacStdDataDir = file.path(gdacDir, sprintf("stddata__%s", gdacStdDataDate))

gdacAnalysesDate = "2014_10_17"
gdacAnalysesDateNoUb = gsub("_", "", gdacAnalysesDate)
gdacAnalysesDir = file.path(gdacDir, sprintf("analyses__%s", gdacAnalysesDate))


##
## DF for combined cancer gene list
##
#cancerGenesFile = file.path(tableDir, "Cancer_Genes.txt")
#cancerGenesGrpModDf = read.delim(cancerGenesFile, header = T, as.is = T)
#cancerGenes = as.character(unique(cancerGenesGrpModDf$gene))


##
## Read TCGA Broad GDAC data
#ncanType = ifelse(canType == "CRC", "COADREAD", canType)
#gdacStdDataCanDir = file.path(gdacStdDataDir, ncanType, gdacStdDataDateNoUb)
#gdacAnalysesCanDir = file.path(gdacAnalysesDir, ncanType, gdacAnalysesDateNoUb)
#canVcfDir = file.path(germDir, "VCF/CANCERS", canType)  

## MutSig SMG list
#mutSigGenesFile = file.path(gdacAnalysesCanDir,
                            #sprintf("gdac.broadinstitute.org_%s-TP.MutSigNozzleReport1.5.Level_4.%s00.0.0", ncanType, gdacAnalysesDateNoUb),
                            #sprintf("%s-TP.sig_genes.txt", ncanType))
#mutSigGenesDf = read.delim(mutSigGenesFile, header = T, as.is = T)
#mutSigGenesDf$nq = gsub("<", "", mutSigGenesDf$q)
#mutSigGenesDf$nq = as.numeric(mutSigGenesDf$nq)
#mutSigGenes = mutSigGenesDf[mutSigGenesDf$nq < 0.05,]$gene


##
## Read gene-level somatic SCNA file
## : The copy number values in the table are in units of (copy number -2),
##   so that no amplification or deletion is 0, genes with amplifications have positive values, and genes with deletions are negative values.
##
#if (ncanType == "SKCM") {
    #scnaFile = file.path(gdacAnalysesCanDir,
                         #sprintf("gdac.broadinstitute.org_%s-TM.CopyNumber_Gistic2.Level_4.%s00.0.0/all_data_by_genes.txt", ncanType, gdacAnalysesDateNoUb))
#} else {
    #scnaFile = file.path(gdacAnalysesCanDir,
                         #sprintf("gdac.broadinstitute.org_%s-TP.CopyNumber_Gistic2.Level_4.%s00.0.0/all_data_by_genes.txt", ncanType, gdacAnalysesDateNoUb))
#}
#scnaDf = read.delim(scnaFile, header = T, as.is = T)
#rownames(scnaDf) = scnaDf$Gene.Symbol
#scnaDf = scnaDf[,4:ncol(scnaDf)]
#scnaDf = scnaDf + 2
#colnames(scnaDf) = sapply(colnames(scnaDf), function(x) { paste(strsplit(x, ".", fixed = T)[[1]][1:3], collapse = "_") })
#scnaPatientIds = colnames(scnaDf)


##
## Read methylation data for promoter regions
##
#methylMergedFile = file.path(tcgaCanDir, sprintf("Methylation-Tumor-%s.txt", canType))
#methylMergedDf = read.delim(methylMergedFile, header = T, as.is = T)
#rownames(methylMergedDf) = methylMergedDf[, "Gene"]


##
## Read gene expression data
##
#geneExpMergedFile = file.path(tcgaCanDir,  sprintf("Gene_Expression-Tumor-%s.txt", canType))
#geneExpMergedDf = read.delim(geneExpMergedFile, header = T, row.names = 1, as.is = T)


##
## Patient IDs falling in the intersetion between data sets
##
#sharedExpScnaMethylPatients = mixedsort(intersect(intersect(colnames(geneExpMergedDf), colnames(scnaDf)), colnames(methylMergedDf)))
#sharedSnpExpScnaMethylEurPatients = mixedsort(intersect(eurSnpPatientIds, sharedExpScnaMethylPatients))


##
## Read VCFs and create genotypes file for prioritization
##
vcfFiles = mixedsort(Sys.glob(file.path(panVcfDir, "*.sanitized.recurrent.vcf.gz")))
vcfFile = vcfFiles[17]

for (vcfFile in vcfFiles) {
    cat(sprintf("Reading %s ...\n", vcfFile))
    vcf = readVcf(vcfFile, "hg19")

    # Recurrent SNVs
    #print(max(info(vcf)$SN))
    print(sort(unique(info(vcf)$SN)))

    regSnvs = vcf[(info(vcf)$SN == 6)]
    info(regSnvs)

    regSnvs = vcf[(info(vcf)$SN == 6) &
                  (elementLengths(info(vcf)$wgEncodeRegDnaseClusteredV3) > 0 |
                   elementLengths(info(vcf)$wgEncodeRegTfbsClusteredWithCellsV3) > 0 |
                   elementLengths(info(vcf)$distalDhsToPromoterDhs) > 0 |
                   elementLengths(info(vcf)$dhsToGeneExpression) > 0 |
                   elementLengths(info(vcf)$fantom5PermissiveEnhancers) > 0 |
                   elementLengths(info(vcf)$fantom5EnhancerTssAssociations) > 0)]

    max(info(regSnvs)$SN)

    nrow(regSnvs)
    info(regSnvs)

    # One-off SNVs
    oneoffSnvs = vcf[info(vcf)$SN == 1]
    nrow(oneoffSnvs)

    #info(novelNoRsVcf)[grepl("transcript_ablation", lapply(info(novelNoRsVcf)$CSQ, function(x) x[1])),]
    #info(novelNoRsVcf)[grepl("splice_donor_variant", lapply(info(novelNoRsVcf)$CSQ, function(x) x[1])),]
    #info(novelNoRsVcf)[grepl("splice_acceptor_variant", lapply(info(novelNoRsVcf)$CSQ, function(x) x[1])),]
    #info(novelNoRsVcf)[grepl("stop_gained", lapply(info(novelNoRsVcf)$CSQ, function(x) x[1])),]
    #info(novelNoRsVcf)[grepl("frameshift_variant", lapply(info(novelNoRsVcf)$CSQ, function(x) x[1])),]
    #info(novelNoRsVcf)[grepl("stop_lost", lapply(info(novelNoRsVcf)$CSQ, function(x) x[1])),]
    #info(novelNoRsVcf)[grepl("initiator_codon_variant", lapply(info(novelNoRsVcf)$CSQ, function(x) x[1])),]
    #info(novelNoRsVcf)[grepl("transcript_amplification", lapply(info(novelNoRsVcf)$CSQ, function(x) x[1])),]
    #info(novelNoRsVcf)[grepl("inframe_insertion", lapply(info(novelNoRsVcf)$CSQ, function(x) x[1])),]
    #info(novelNoRsVcf)[grepl("inframe_deletion", lapply(info(novelNoRsVcf)$CSQ, function(x) x[1])),]
    #info(novelNoRsVcf)[grepl("missense_variant", lapply(info(novelNoRsVcf)$CSQ, function(x) x[1])),]
    #info(novelNoRsVcf)[grepl("splice_region_variant", lapply(info(novelNoRsVcf)$CSQ, function(x) x[1])),]
    #info(novelNoRsVcf)[grepl("incomplete_terminal_codon_variant", lapply(info(novelNoRsVcf)$CSQ, function(x) x[1])),]
    #info(novelNoRsVcf)[grepl("stop_retained_variant", lapply(info(novelNoRsVcf)$CSQ, function(x) x[1])),]
    #info(novelNoRsVcf)[grepl("synonymous_variant", lapply(info(novelNoRsVcf)$CSQ, function(x) x[1])),]
    #info(novelNoRsVcf)[grepl("coding_sequence_variant", lapply(info(novelNoRsVcf)$CSQ, function(x) x[1])),]
    #info(novelNoRsVcf)[grepl("TFBS_ablation", lapply(info(novelNoRsVcf)$CSQ, function(x) x[1])),]
    #info(novelNoRsVcf)[grepl("TFBS_amplification", lapply(info(novelNoRsVcf)$CSQ, function(x) x[1])),]
    #info(novelNoRsVcf)[grepl("TF_binding_site_variant", lapply(info(novelNoRsVcf)$CSQ, function(x) x[1])),]
    #info(novelNoRsVcf)[grepl("regulatory_region_variant", lapply(info(novelNoRsVcf)$CSQ, function(x) x[1])),]
    #info(novelNoRsVcf)[grepl("intergenic_variant", lapply(info(novelNoRsVcf)$CSQ, function(x) x[1])),]
}



    ##
    ## One to one assocation tests between genotype group and expression change sof nearby/associated gene(s)
    ##
    upstreamCgSnpsSigExpDf = data.frame()

    for (i in 1:nrow(upstreamCgSnpsDf)) {
        print(i)
        r = upstreamCgSnpsDf[i,]

        mutSampleIds = colnames(r)[grepl("'\\S{1}/\\S{1}'", r, perl = TRUE) & 
                                   !grepl("0/0", r, perl = TRUE) & 
                                   !grepl("\\./\\.", r, perl = TRUE)]
        mutPatientIds = gsub("(TCGA_\\S{2}_\\S{4})_\\S{3}", "\\1", mutSampleIds)
        wdtSampleIds = colnames(r)[grepl("0/0", r, perl = TRUE)]
        wdtPatientIds = gsub("(TCGA_\\S{2}_\\S{4})_\\S{3}", "\\1", wdtSampleIds)

        nearbyGeneNames = strsplit(gsub("\\(.*?\\)", "", r$feature), ",")[[1]]
        corrDhsDhsGeneNames = strsplit(gsub("\\(.*?\\)", "", gsub("Name=", "", r$distal_dhs_to_promoter_dhs)), ",")[[1]]
        corrDhsRnaGeneNames = strsplit(gsub("\\(.*?\\)", "", gsub("Name=", "", r$dhs_to_gene_expression)), ",")[[1]]
        geneNames = unique(na.omit(c(nearbyGeneNames, corrDhsDhsGeneNames, corrDhsRnaGeneNames)))

        for (geneName in geneNames) {
            if (geneName == "CAAP1") geneName = "C9orf82"
            varStem = sprintf("%s_%d_%s_%s_%s", r$chrom, r$pos, r$ref, r$alt, geneName)
            figStem = sprintf("chr%s:%d:%s>%s", r$chrom, r$pos, r$ref, r$alt)
            selGeneExps = geneExpMergedDf[grep(paste(geneName, "|", sep=""), rownames(geneExpMergedDf), fixed = T),]
            if (nrow(selGeneExps) != 1) {
                print(sprintf("No expression data for %s", geneName))
                next
            }
      
            mutGeneExps = selGeneExps[, which(gsub("(TCGA_\\S{2}_\\S{4})_\\S{3}", "\\1", colnames(selGeneExps)) %in% mutPatientIds)]
            wdtGeneExps = selGeneExps[, which(gsub("(TCGA_\\S{2}_\\S{4})_\\S{3}", "\\1", colnames(selGeneExps)) %in% wdtPatientIds)]
            
            if (length(mutGeneExps) == 0) {
                print(sprintf("No %s mutant samples in expression data", geneName))
                next
            }

            if (length(wdtGeneExps) == 0) {
                print(sprintf("No %s wildtype samples in expression data", geneName))
                next
            }

            ## RSEM
            selGeneExpMutVsWdtDf = data.frame()
            selGeneExpMutVsWdtDf = rbind(selGeneExpMutVsWdtDf, data.frame(type=rep(figStem, length(mutGeneExps)), rsem = as.numeric(mutGeneExps)))
            selGeneExpMutVsWdtDf = rbind(selGeneExpMutVsWdtDf, data.frame(type=rep("WT", length(wdtGeneExps)), rsem = as.numeric(wdtGeneExps)))

            wtPvalueGeneExpMutVsWdt = my.wilcox.test.p.value(as.numeric(mutGeneExps), as.numeric(wdtGeneExps))
            ttPvalueGeneExpMutVsWdt = my.t.test.p.value(as.numeric(mutGeneExps), as.numeric(wdtGeneExps))

            pvalueCondGeneExpMutVsWdt = (wtPvalueGeneExpMutVsWdt < 0.05 | ttPvalueGeneExpMutVsWdt < 0.05)

            if (!is.na(pvalueCondGeneExpMutVsWdt) & pvalueCondGeneExpMutVsWdt == TRUE) {
                mutVsWdtExpPlot = ggplot(selGeneExpMutVsWdtDf, aes(factor(type), rsem)) +
                ggtitle(sprintf("%s\n", geneName)) +
                geom_boxplot() + geom_jitter() +
                theme(plot.title   = element_text(size = baseFontSize + 2, face="bold"),
                      axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
                      axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
                      axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 0, hjust = NULL),
                      axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(color = "black", fill="white"),
                      legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
                      legend.text  = element_text(size = baseFontSize, family="sans"),
                      legend.direction = "horizontal",
                      legend.position = "bottom") +
                                   scale_x_discrete(name=sprintf("\nWilcoxon signed-rank test: %f\nStudent's t-test: %f", wtPvalueGeneExpMutVsWdt, ttPvalueGeneExpMutVsWdt)) +
                                   scale_y_continuous(name = "RSEM\n")
                                   mutVsWdtExpPlotFile = file.path(expFigDir, sprintf("%s-MutVsWdt-RSEM.pdf", varStem))
                                   ggsave(filename = mutVsWdtExpPlotFile, plot = mutVsWdtExpPlot, width = 4, height = 5)
            }

            #
            selGeneExpMutVsTotDf = data.frame()
            selGeneExpMutVsTotDf = rbind(selGeneExpMutVsTotDf, data.frame(type=rep(figStem, length(mutGeneExps)), log2FoldChange = as.numeric(mutGeneExps)))
            selGeneExpMutVsTotDf = rbind(selGeneExpMutVsTotDf, data.frame(type=rep("Total", length(selGeneExps)), log2FoldChange = as.numeric(selGeneExps)))

            wtPvalueGeneExpMutVsTot = my.wilcox.test.p.value(as.numeric(mutGeneExps), as.numeric(selGeneExps))
            ttPvalueGeneExpMutVsTot = my.t.test.p.value(as.numeric(mutGeneExps), as.numeric(selGeneExps))

            pvalueCondGeneExpMutVsTot = (wtPvalueGeneExpMutVsTot < 0.05 | ttPvalueGeneExpMutVsTot < 0.05)

            if (!is.na(pvalueCondGeneExpMutVsTot) & pvalueCondGeneExpMutVsTot == TRUE) {
                mutVsTotExpPlot = ggplot(selGeneExpMutVsTotDf, aes(factor(type), log2FoldChange)) +
                ggtitle(sprintf("%s\n", geneName)) +
                geom_boxplot() + geom_jitter() +
                theme(plot.title   = element_text(size = baseFontSize + 2, face="bold"),
                      axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
                      axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
                      axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 0, hjust = NULL),
                      axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(color = "black", fill="white"),
                      legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
                      legend.text  = element_text(size = baseFontSize, family="sans"),
                      legend.direction = "horizontal",
                      legend.position = "bottom") +
                                   scale_x_discrete(name=sprintf("\nWilcoxon signed-rank test: %f\nStudent's t-test: %f", wtPvalueGeneExpMutVsTot, ttPvalueGeneExpMutVsTot)) +
                                   scale_y_continuous(name = "RSEM\n")
                                   mutVsTotExpPlotFile = file.path(expFigDir, sprintf("%s-MutVsTot-RSEM.pdf", varStem))
                                   ggsave(filename = mutVsTotExpPlotFile, plot = mutVsTotExpPlot, width = 4, height = 5)
            }

            ## Log2 Fold Change of RSEM
            selGeneExpLog2Fcs = geneExpLog2FcPrimDf[grep(paste(geneName, "|", sep=""), rownames(geneExpLog2FcPrimDf), fixed = T),]
            mutGeneExpLog2Fcs = selGeneExpLog2Fcs[, which(gsub("(TCGA_\\S{2}_\\S{4})_\\S{3}", "\\1", colnames(selGeneExpLog2Fcs)) %in% mutPatientIds)]
            wdtGeneExpLog2Fcs = selGeneExpLog2Fcs[, which(gsub("(TCGA_\\S{2}_\\S{4})_\\S{3}", "\\1", colnames(selGeneExpLog2Fcs)) %in% wdtPatientIds)]

            selGeneExpLog2FcMutVsWdtDf = data.frame()
            selGeneExpLog2FcMutVsWdtDf = rbind(selGeneExpLog2FcMutVsWdtDf, 
                                               data.frame(type=rep(figStem, length(mutGeneExpLog2Fcs)), log2FoldChange = as.numeric(mutGeneExpLog2Fcs)))
            selGeneExpLog2FcMutVsWdtDf = rbind(selGeneExpLog2FcMutVsWdtDf, 
                                               data.frame(type=rep("WT", length(wdtGeneExpLog2Fcs)), log2FoldChange = as.numeric(wdtGeneExpLog2Fcs)))

            wtPvalueGeneExpLog2FcMutVsWdt = my.wilcox.test.p.value(as.numeric(mutGeneExpLog2Fcs), as.numeric(wdtGeneExpLog2Fcs))
            ttPvalueGeneExpLog2FcMutVsWdt = my.t.test.p.value(as.numeric(mutGeneExpLog2Fcs), as.numeric(wdtGeneExpLog2Fcs))

            pvalueCondGeneExpLog2FcMutVsWdt = (wtPvalueGeneExpLog2FcMutVsWdt < 0.05 | ttPvalueGeneExpLog2FcMutVsWdt < 0.05)

            if (!is.na(pvalueCondGeneExpLog2FcMutVsWdt) & pvalueCondGeneExpLog2FcMutVsWdt == TRUE) {
                mutVsWdtExpLog2FcPlot = ggplot(selGeneExpLog2FcMutVsWdtDf, aes(factor(type), log2FoldChange)) +
                ggtitle(sprintf("%s\n", geneName)) +
                geom_boxplot() + geom_jitter() +
                theme(plot.title   = element_text(size = baseFontSize + 2, face="bold"),
                      axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
                      axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
                      axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 0, hjust = NULL),
                      axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(color = "black", fill="white"),
                      legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
                      legend.text  = element_text(size = baseFontSize, family="sans"),
                      legend.direction = "horizontal",
                      legend.position = "bottom") +
                                   scale_x_discrete(name=sprintf("\nWilcoxon signed-rank test: %f\nStudent's t-test: %f", wtPvalueGeneExpLog2FcMutVsWdt, ttPvalueGeneExpLog2FcMutVsWdt)) +
                                   scale_y_continuous(name = "Log2 fold change of RSEM\n")
                                   mutVsWdtExpLog2FcPlotFile = file.path(expFigDir, sprintf("%s-MutVsWdt-Log2FoldChangeRSEM.pdf", varStem))
                                   ggsave(filename = mutVsWdtExpLog2FcPlotFile, plot = mutVsWdtExpLog2FcPlot, width = 4, height = 5)
            }

            #
            selGeneExpLog2FcMutVsTotDf = data.frame()
            selGeneExpLog2FcMutVsTotDf = rbind(selGeneExpLog2FcMutVsTotDf, data.frame(type=rep(figStem, length(mutGeneExpLog2Fcs)), log2FoldChange = as.numeric(mutGeneExpLog2Fcs)))
            selGeneExpLog2FcMutVsTotDf = rbind(selGeneExpLog2FcMutVsTotDf, data.frame(type=rep("Total", length(selGeneExpLog2Fcs)), log2FoldChange = as.numeric(selGeneExpLog2Fcs)))

            wtPvalueGeneExpLog2FcMutVsTot = my.wilcox.test.p.value(as.numeric(mutGeneExpLog2Fcs), as.numeric(selGeneExpLog2Fcs))
            ttPvalueGeneExpLog2FcMutVsTot = my.t.test.p.value(as.numeric(mutGeneExpLog2Fcs), as.numeric(selGeneExpLog2Fcs))

            pvalueCondGeneExpLog2FcMutVsTot = (wtPvalueGeneExpLog2FcMutVsTot < 0.05 | ttPvalueGeneExpLog2FcMutVsTot < 0.05)

            if (!is.na(pvalueCondGeneExpLog2FcMutVsTot) & pvalueCondGeneExpLog2FcMutVsTot == TRUE) {
                mutVsTotExpLog2FcPlot = ggplot(selGeneExpLog2FcMutVsTotDf, aes(factor(type), log2FoldChange)) +
                ggtitle(sprintf("%s\n", geneName)) +
                geom_boxplot() + geom_jitter() +
                theme(plot.title   = element_text(size = baseFontSize + 2, face="bold"),
                      axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
                      axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
                      axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 0, hjust = NULL),
                      axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(color = "black", fill="white"),
                      legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
                      legend.text  = element_text(size = baseFontSize, family="sans"),
                      legend.direction = "horizontal",
                      legend.position = "bottom") +
                                   scale_x_discrete(name=sprintf("\nWilcoxon signed-rank test: %f\nStudent's t-test: %f", 
                                                                 wtPvalueGeneExpLog2FcMutVsTot, 
                                                                 ttPvalueGeneExpLog2FcMutVsTot)) +
                                   scale_y_continuous(name = "Log2 fold change of RSEM\n")
                                   mutVsTotExpLog2FcPlotFile = file.path(expFigDir, sprintf("%s-MutVsTot-Log2FoldChangeRSEM.pdf", varStem))
                                   ggsave(filename = mutVsTotExpLog2FcPlotFile, plot = mutVsTotExpLog2FcPlot, width = 4, height = 5)
            }

            if (!is.na(pvalueCondGeneExpMutVsWdt) & pvalueCondGeneExpMutVsWdt == TRUE) {
                upstreamCgSnpsSigExpDf = rbind(upstreamCgSnpsSigExpDf,
                                                 data.frame(r, 
                                                            gene = geneName, 
                                                            mean_mutant_gene_expression = mean(as.numeric(mutGeneExps)),
                                                            mean_wildtype_gene_expression = mean(as.numeric(wdtGeneExps)),
                                                            mean_total_gene_expression = mean(as.numeric(selGeneExps)),
                                                            median_mutant_gene_expression = median(as.numeric(mutGeneExps)),
                                                            median_wildtype_gene_expression = median(as.numeric(wdtGeneExps)),
                                                            median_total_gene_expression = median(as.numeric(selGeneExps)),
                                                            mean_mutant_gene_expression_log2_fold_change = mean(as.numeric(mutGeneExpLog2Fcs)),
                                                            mean_wildtype_gene_expression_log2_fold_change = mean(as.numeric(wdtGeneExpLog2Fcs)),
                                                            mean_total_gene_expression_log2_fold_change = mean(as.numeric(selGeneExpLog2Fcs)),
                                                            median_mutant_gene_expression_log2_fold_change = median(as.numeric(mutGeneExpLog2Fcs)),
                                                            median_wildtype_gene_expression_log2_fold_change = median(as.numeric(wdtGeneExpLog2Fcs)),
                                                            median_total_gene_expression_log2_fold_change = median(as.numeric(selGeneExpLog2Fcs)),
                                                            wt_pvalue_gene_expression_mutant_vs_wildtype = wtPvalueGeneExpMutVsWdt,
                                                            tt_pvalue_gene_expression_mutant_vs_wildtype = ttPvalueGeneExpMutVsWdt,
                                                            wt_pvalue_gene_expression_mutant_vs_total = wtPvalueGeneExpMutVsTot,
                                                            tt_pvalue_gene_expression_mutant_vs_total = ttPvalueGeneExpMutVsTot,
                                                            wt_pvalue_gene_expression_log2_fold_change_mutant_vs_wildtype = wtPvalueGeneExpLog2FcMutVsWdt,
                                                            tt_pvalue_gene_expression_log2_fold_change_mutant_vs_wildtype = ttPvalueGeneExpLog2FcMutVsWdt,
                                                            wt_pvalue_gene_expression_log2_fold_change_mutant_vs_total = wtPvalueGeneExpLog2FcMutVsTot,
                                                            tt_pvalue_gene_expression_log2_fold_change_mutant_vs_total = ttPvalueGeneExpLog2FcMutVsTot))
            }
        }
    }

    upstreamCgSnpsSigExpDf2 = subset(upstreamCgSnpsSigExpDf,
                                       (!is.na(wt_pvalue_gene_expression_mutant_vs_wildtype) & wt_pvalue_gene_expression_mutant_vs_wildtype <= 0.01) &
                                       (!is.na(tt_pvalue_gene_expression_mutant_vs_wildtype) & tt_pvalue_gene_expression_mutant_vs_wildtype <= 0.01) &
                                       tg_all_af1 < 0.05)

    upstreamCgSnpsSigExpDf2 = upstreamCgSnpsSigExpDf2[order(-upstreamCgSnpsSigExpDf2$all_af1),]
    #upstreamCgSnpsSigExpDf2 = upstreamCgSnpsSigExpDf2[order(upstreamCgSnpsSigExpDf2$wt_pvalue_gene_expression_mutant_vs_wildtype),]

    #upstreamCgSnpsSigExpDf2 = subset(upstreamCgSnpsSigExpDf,
                                       #(!is.na(wt_pvalue_gene_expression_mutant_vs_wildtype) & wt_pvalue_gene_expression_mutant_vs_wildtype <= 0.05) |
                                       #(!is.na(tt_pvalue_gene_expression_mutant_vs_wildtype) & tt_pvalue_gene_expression_mutant_vs_wildtype <= 0.05))

    upstreamCgSnpsSigExpDf2[c(31:40),]

    nrow(upstreamCgSnpsSigExpDf2)
    head(upstreamCgSnpsSigExpDf2)






    # Promoter DHS
    novel05SnpsPromoDf = subset(novel05SnpsSrtDf, !is.na(promoter_dhs_master_known) & region != "upstream")
    nrow(novel05SnpsPromoDf)
    head(novel05SnpsPromoDf)
    novel05SnpsPromoterDHSFile = file.path(tableDir,  sprintf("GCC-Germline-%s-PromoterDHS.txt", canType))
    write.table(novel05SnpsPromoterDHSDf, file=novel05SnpsPromoterDHSFile, sep="\t", row.names=F, col.names=T, quote=FALSE)

    # Distal DHS (DHS-DHS)
    novel05SnpsDhsDhsDf = subset(novel05SnpsDf, !is.na(distal_dhs_to_promoter_dhs))
    nrow(novel05SnpsDhsDhsDf)
    head(novel05SnpsDf,2)

    # Distal DHS (DHS-mRNA)
    novel05SnpsDhsExpDf = subset(novel05SnpsDf, !is.na(dhs_to_gene_expression))


    # Upstream or downstream

    ## All known SNPs
    knownSnpsDf = sqldf(sprintf("SELECT *
                                FROM %s.snp_ext_annotations 
                                WHERE (chrom != 'X' AND chrom != 'Y') AND
                                (id != '.' OR rsid1 != '' OR dbsnp137 IS NOT NULL OR dbsnp138 IS NOT NULL)", canType), dbname = sqlite3Db)

    # Known common SNPs with MAF >= 0.05 (0.05 <= AF <= 0.95)
    knownCommonSnpsDf = sqldf(sprintf("SELECT *
                                        FROM %s.snp_annotations 
                                        WHERE  (chrom != 'X' AND chrom != 'Y') AND
                                        (id != '.' OR rsid1 != '' OR dbsnp137 IS NOT NULL OR dbsnp138 IS NOT NULL) AND 
                                        NOT ((tg_all_af1 < 0.05 AND tg_asn_af1 < 0.05 AND tg_amr_af1 < 0.05 AND tg_afr_af1 < 0.05 AND tg_eur_af1 < 0.05) OR
                                        (tg_all_af1 > 0.95 AND tg_asn_af1 > 0.95 AND tg_amr_af1 > 0.95 AND tg_afr_af1 > 0.95 AND tg_eur_af1 > 0.95))", canType), dbname = sqlite3Db)

    # Known rare SNPs with MAF < 0.05 (AF1 < 0.05 or AF > 0.95)
    #knownRareSnpsDf = sqldf(sprintf("SELECT *
                                        #FROM %s.snp_annotations 
                                        #WHERE  (chrom != 'X' AND chrom != 'Y') AND
                                        #(id != '.' OR rsid1 != '' OR dbsnp137 IS NOT NULL OR dbsnp138 IS NOT NULL) AND 
                                        #(tg_all_af1 IS NOT NULL AND tg_asn_af1 IS NOT NULL AND tg_amr_af1 IS NOT NULL AND tg_afr_af1 IS NOT NULL AND tg_eur_af1 IS NOT NULL) AND
                                        #(tg_all_af1 != '' AND tg_asn_af1 != '' AND tg_amr_af1 != '' AND tg_afr_af1 != '' AND tg_eur_af1 != '') AND
                                        #(tg_all_af1 > 0 AND tg_asn_af1 >= 0 AND tg_amr_af1 >= 0 AND tg_afr_af1 >= 0 AND tg_eur_af1 >= 0) AND
                                        #((tg_all_af1 < 0.05 AND tg_asn_af1 < 0.05 AND tg_amr_af1 < 0.05 AND tg_afr_af1 < 0.05 AND tg_eur_af1 < 0.05) OR
                                        #(tg_all_af1 > 0.95 AND tg_asn_af1 > 0.95 AND tg_amr_af1 > 0.95 AND tg_afr_af1 > 0.95 AND tg_eur_af1 > 0.95))", canType), dbname = sqlite3Db)

    knownRareSnpsDf = sqldf(sprintf("SELECT *
                                    FROM %s.snp_annotations 
                                    WHERE (chrom != 'X' AND chrom != 'Y') AND
                                    (tg_all_af1 IS NOT NULL AND  tg_all_af1 != '' AND tg_all_af1 > 0) AND 
                                    (tg_all_af1 < 0.05 OR tg_all_af1 > 0.95)", canType), dbname = sqlite3Db)

    head(knownRareSnpsDf)
    nrow(knownRareSnpsDf)
    summary(knownRareSnpsDf$cadd_phred1)
    summary(knownRareSnpsDf$cadd_raw1)

    # Calculate adjusted p-values and create a separate table for known rare SNPs
    tableColNames = sqldf("pragma table_info(snp_annotations)", dbname = sqlite3Db)$name
    tablePvalueCols = grep("adj", grep("pval", tableColNames, value=T), invert=T, value=T)
    for (tablePvalueCol in tablePvalueCols) {
        adjPvalueCol = gsub("pval", "adj_pval", tablePvalueCol)
        knownRareSnpsDf[, adjPvalueCol] = p.adjust(knownRareSnpsDf[, tablePvalueCol], "fdr")
    }

    sqldf(sprintf("DROP TABLE IF EXISTS %s.known_rare_snps; CREATE TABLE %s.known_rare_snps AS SELECT * FROM main.knownRareSnpsDf ", canType, canType), dbname = sqlite3Db)


    orKnownRareSnpsDf = sqldf(sprintf("SELECT *
                                      FROM main.knownRareSnpsDf
                                      WHERE (alt2 = '' AND alt3 = '') AND (fisher_adj_pval_all_ac1 < 0.01 AND chisq_adj_pval_all_ac1 < 0.01) AND
                                            ((tg_all_af1 > 0.95 AND all_af1 < 0.95) OR (tg_all_af1 < 0.05 AND all_af1 > 0.05))", canType), dbname = sqlite3Db)
    
    orKnownRareSnpsDf = orKnownRareSnpsDf[order(-orKnownRareSnpsDf$all_af1),]
    summary(orKnownRareSnpsDf$cadd_phred1)
    summary(orKnownRareSnpsDf$cadd_raw1)

    head(orKnownRareSnpsDf)
    tail(orKnownRareSnpsDf)
    nrow(orKnownRareSnpsDf)

    ## ChromHMMs for SNPs
    #knownSnpsBroadHmmDf = data.frame()
    #broadHmmCols = grep("broad_hmm", colnames(knownSnpsDf), value=T)
    #for (broadHmmCol in broadHmmCols) {
        #hmmTmpGrpDf = sqldf(sprintf("SELECT %s, COUNT(*) AS mutCnt FROM main.knownSnpsDf GROUP BY %s", broadHmmCol, broadHmmCol))
        #hmmTmpGrpDf[,broadHmmCol] = gsub("Name=", "", hmmTmpGrpDf[,broadHmmCol], fixed=T)
        #hmmTmpCellType = toupper(strsplit(broadHmmCol, "_", fixed=T)[[1]][5])
        #hmmTmpGrpDf$cellType = hmmTmpCellType
        #colnames(hmmTmpGrpDf)[1] = "state"
        #hmmTmpGrpDf2 = sqldf('SELECT a.*, h.size, h.color
                                    #FROM main.hmmTmpGrpDf AS a 
                                    #LEFT JOIN main.chromHmmBedDf AS h 
                                    #ON a.state = h.state AND a.cellType = h.cellType')

        #hmmTmpGrpDf2$mutRate = hmmTmpGrpDf2$mutCnt / hmmTmpGrpDf2$size
        #hmmTmpGrpDf2$state = gsub("_", " ", hmmTmpGrpDf2$state)
        #hmmTmpGrpDf2$state = factor(hmmTmpGrpDf2$state, levels=ord)
        #hmmTmpGrpDf2 = hmmTmpGrpDf2[order(hmmTmpGrpDf2$state),]

        #hmmTmpGrpDf2$cumSumMutRate = cumsum(hmmTmpGrpDf2$mutRate)
        #hmmTmpGrpDf2$cumSumMutRateMidPoint = NA
        #for (i in 1:nrow(hmmTmpGrpDf2)) {
            #if (i == 1) {
                #hmmTmpGrpDf2[i,]$cumSumMutRateMidPoint = hmmTmpGrpDf2[i,]$cumSumMutRate / 2
            #} else {
                #hmmTmpGrpDf2[i,]$cumSumMutRateMidPoint = hmmTmpGrpDf2[i-1,]$cumSumMutRate + ((hmmTmpGrpDf2[i,]$cumSumMutRate - hmmTmpGrpDf2[i-1,]$cumSumMutRate) / 2)
            #}
        #}

        #hmmTmpGrpDf2$cumSumMutCnt = cumsum(hmmTmpGrpDf2$mutCnt)
        #hmmTmpGrpDf2$cumSumMutCntMidPoint = NA
        #for (i in 1:nrow(hmmTmpGrpDf2)) {
            #if (i == 1) {
                #hmmTmpGrpDf2[i,]$cumSumMutCntMidPoint = hmmTmpGrpDf2[i,]$cumSumMutCnt / 2
            #} else {
                #hmmTmpGrpDf2[i,]$cumSumMutCntMidPoint = hmmTmpGrpDf2[i-1,]$cumSumMutCnt + ((hmmTmpGrpDf2[i,]$cumSumMutCnt - hmmTmpGrpDf2[i-1,]$cumSumMutCnt) / 2)
            #}
        #}

        #hmmTmpGrpDf2 = hmmTmpGrpDf2[!is.na(hmmTmpGrpDf2$state),]
        #totMutCnt = sum(as.numeric(hmmTmpGrpDf2$mutCnt))
        #totStateSize = sum(as.numeric(hmmTmpGrpDf2$size))
        #avgMutCnt = totMutCnt / totStateSize

        #for (i in 1:nrow(hmmTmpGrpDf2)) {
            #expMutCnt = hmmTmpGrpDf2[i,]$size * avgMutCnt
            #poisPvalue = poisson.test(hmmTmpGrpDf2[i,]$mutCnt, expMutCnt)$p.value
            #hmmTmpGrpDf2[i, "poisPvalue"] = poisPvalue
            #ct = matrix(c(hmmTmpGrpDf2[i,]$mutCnt, as.numeric(hmmTmpGrpDf2[i,]$size) - hmmTmpGrpDf2[i,]$mutCnt,
                          #sum(hmmTmpGrpDf2[-i,]$mutCnt), sum(as.numeric(hmmTmpGrpDf2[-i,]$size)) - sum(hmmTmpGrpDf2[-i,]$mutCnt)),
                        #2,2, dimnames=list(mutStatus=c("mut", "nonMut"),sites=c("A","nonA")))
            ##c = chisq.test(ct)
            ##hmmTmpGrpDf2[i, "chisqPvalue"] = c$p.value
            #g = g.test(ct)
            #hmmTmpGrpDf2[i, "gPvalue"] = g$p.value
        #}
        #hmmTmpGrpDf2$poisAdjPvalue= p.adjust(hmmTmpGrpDf2$poisPvalue, method="fdr")
    #}


    # Mutation spectra of known SNPs
    knownSnpsDf2 = knownSnpsDf
    knownSnpsDf2$new_ref = knownSnpsDf2$ref
    knownSnpsDf2$new_alt = knownSnpsDf2$alt

    knownSnpsDf2[knownSnpsDf2$ref == "G", "new_ref"] = "C"
    knownSnpsDf2[knownSnpsDf2$ref == "G" & knownSnpsDf2$alt == "A", "new_alt"] = "T"
    knownSnpsDf2[knownSnpsDf2$ref == "G" & knownSnpsDf2$alt == "C", "new_alt"] = "G"
    knownSnpsDf2[knownSnpsDf2$ref == "G" & knownSnpsDf2$alt == "T", "new_alt"] = "A"

    knownSnpsDf2[knownSnpsDf2$ref == "A", "new_ref"] = "T"
    knownSnpsDf2[knownSnpsDf2$ref == "A" & knownSnpsDf2$alt == "G", "new_alt"] = "C"
    knownSnpsDf2[knownSnpsDf2$ref == "A" & knownSnpsDf2$alt == "C", "new_alt"] = "G"
    knownSnpsDf2[knownSnpsDf2$ref == "A" & knownSnpsDf2$alt == "T", "new_alt"] = "A"

    knownSnpsDf2$Substitution_Class = sprintf("%s>%s", knownSnpsDf2$new_ref, knownSnpsDf2$new_alt)
    knownSnpsDf2$Substitution_Class = factor(knownSnpsDf2$Substitution_Class, level = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))

    knownSnpsSubGrpDf = sqldf('SELECT Substitution_Class, COUNT(*) AS Count FROM knownSnpsDf2 GROUP BY Substitution_Class')
    knownSnpsSubGrpDf$Type = "Known SNPs"

    # PhastCons of known SNPs
    knownSnpsDf2$phast_cons_elements46way_score = as.numeric(gsub(".*lod=(\\d+).*", "\\1", knownSnpsDf2$phast_cons_elements46way))
    knownSnpsDf2$phast_cons_elements46way_placental_score = as.numeric(gsub(".*lod=(\\d+).*", "\\1", knownSnpsDf2$phast_cons_elements46way_placental))
    knownSnpsDf2$phast_cons_elements46way_primates_score = as.numeric(gsub(".*lod=(\\d+).*", "\\1", knownSnpsDf2$phast_cons_elements46way_primates))
    boxplot(knownSnpsDf2$phast_cons_elements46way_primates_score)
    #head(knownSnpsDf2$phast_cons_elements46way, 100)
    summary(knownSnpsDf2$phast_cons_elements46way_score)

    # CADD score of known SNPs
    head(knownSnpsDf2$cadd_raw1)
    head(knownSnpsDf2$cadd_phred1)

    # Over-represented known SNPs
    orSnpsDf = sqldf(sprintf("SELECT *
                            FROM %s.snp_ext_annotations 
                            WHERE (chrom != 'X' AND chrom != 'Y') AND
                                (id != '.' OR rsid1 != '' OR dbsnp137 IS NOT NULL OR dbsnp138 IS NOT NULL) AND 
                                (all_af1 > tg_all_af1 AND all_af1 >= 0.05 AND tg_all_af1 < 0.05 AND tg_eur_af1 < 0.05) AND
                                (fisher_adj_pval_all_ac1 < 0.01 AND fisher_adj_pval_all_ac1 != '') AND
                                (chisq_adj_pval_all_ac1 < 0.01 AND chisq_adj_pval_all_ac1 != '')", canType), dbname = sqlite3Db)

    nrow(orSnpsDf)
    summary(orSnpsDf$all_af1)

    # Regional breakdown of over-represented known SNPs
    orSnpsRegGrpDf = sqldf("SELECT region, COUNT(*) AS cnt FROM orSnpsDf GROUP BY region")

    # Mutation spectra of over-represented known SNPs
    orSnpsDf2 = orSnpsDf
    orSnpsDf2$new_ref = orSnpsDf2$ref
    orSnpsDf2$new_alt = orSnpsDf2$alt

    orSnpsDf2[orSnpsDf2$ref == "G", "new_ref"] = "C"
    orSnpsDf2[orSnpsDf2$ref == "G" & orSnpsDf2$alt == "A", "new_alt"] = "T"
    orSnpsDf2[orSnpsDf2$ref == "G" & orSnpsDf2$alt == "C", "new_alt"] = "G"
    orSnpsDf2[orSnpsDf2$ref == "G" & orSnpsDf2$alt == "T", "new_alt"] = "A"

    orSnpsDf2[orSnpsDf2$ref == "A", "new_ref"] = "T"
    orSnpsDf2[orSnpsDf2$ref == "A" & orSnpsDf2$alt == "G", "new_alt"] = "C"
    orSnpsDf2[orSnpsDf2$ref == "A" & orSnpsDf2$alt == "C", "new_alt"] = "G"
    orSnpsDf2[orSnpsDf2$ref == "A" & orSnpsDf2$alt == "T", "new_alt"] = "A"

    orSnpsDf2$Substitution_Class = sprintf("%s>%s", orSnpsDf2$new_ref, orSnpsDf2$new_alt)
    orSnpsDf2$Substitution_Class = factor(orSnpsDf2$Substitution_Class, level = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))

    orSnpsSubGrpDf = sqldf('SELECT Substitution_Class, COUNT(*) AS Count FROM orSnpsDf2 GROUP BY Substitution_Class')
    orSnpsSubGrpDf$Type = "Over-represented known SNPs"

    # PhastCons of over-represented known SNPs
    orSnpsDf2$phast_cons_elements46way_score = as.numeric(gsub(".*lod=(\\d+).*", "\\1", orSnpsDf2$phast_cons_elements46way))
    orSnpsDf2$phast_cons_elements46way_placental_score = as.numeric(gsub(".*lod=(\\d+).*", "\\1", orSnpsDf2$phast_cons_elements46way_placental))
    orSnpsDf2$phast_cons_elements46way_primates_score = as.numeric(gsub(".*lod=(\\d+).*", "\\1", orSnpsDf2$phast_cons_elements46way_primates))
    boxplot(orSnpsDf2$phast_cons_elements46way_score)
    head(orSnpsDf2$phast_cons_elements46way_score)
    summary(orSnpsDf2$phast_cons_elements46way_score)

    # Novel SNPs with AF > 0.1
    novel1SnpsDf = sqldf(sprintf("SELECT *
                                FROM snp_ext_annotations 
                                WHERE (chrom != 'X' AND chrom != 'Y') AND
                                (id == '.' AND rsid1 == '' AND tg_all_af1 == '' AND all_af1 >= 0.1 AND dbsnp137 IS NULL AND dbsnp138 IS NULL)", canType), dbname = sqlite3Db)

    # Novel SNPs with AF < 0.1 and AF >= 0.05
    novel105SnpsDf = sqldf(sprintf("SELECT *
                                FROM snp_ext_annotations 
                                WHERE (chrom != 'X' AND chrom != 'Y') AND
                                (id == '.' AND rsid1 == '' AND tg_all_af1 == '' AND all_af1 < 0.1 AND all_af1 >= 0.05 AND dbsnp137 IS NULL AND dbsnp138 IS NULL)", canType), dbname = sqlite3Db)

    # Novel SNPs with AF < 0.05 and AF >= 0.01
    novel0501SnpsDf = sqldf(sprintf("SELECT *
                                   FROM snp_annotations 
                                   WHERE (chrom != 'X' AND chrom != 'Y') AND
                                   (id == '.' AND rsid1 == '' AND tg_all_af1 == '' AND all_af1 < 0.05 AND all_af1 >= 0.01 AND dbsnp137 IS NULL AND dbsnp138 IS NULL)", canType), dbname = sqlite3Db)

    # Novel SNPs with AF < 0.01
    novel01SnpsDf = sqldf(sprintf("SELECT *
                                   FROM snp_annotations 
                                   WHERE (chrom != 'X' AND chrom != 'Y') AND
                                   (id == '.' AND rsid1 == '' AND tg_all_af1 == '' AND all_af1 < 0.01 AND dbsnp137 IS NULL AND dbsnp138 IS NULL)", canType), dbname = sqlite3Db)

    # Novel SNPs with AF >= 0.05
    #novel05SnpsDf = sqldf(sprintf("SELECT *
                                   #FROM snp_annotations 
                                   #WHERE (chrom != 'X' AND chrom != 'Y') AND
                                   #(id == '.' AND rsid1 == '' AND tg_all_af1 == '' AND all_af1 >= 0.05 AND dbsnp137 IS NULL AND dbsnp138 IS NULL)", canType), dbname = sqlite3Db)

    novel05SnpsDf = sqldf(sprintf("SELECT *
                                   FROM snp_annotations 
                                   WHERE (chrom != 'X' AND chrom != 'Y') AND
                                   (id == '.' AND rsid1 == '' AND tg_all_af1 == '' AND all_af1 >= 0.05 AND dbsnp137 IS NULL AND dbsnp138 IS NULL)", canType), dbname = sqlite3Db)

    novel05SnpsDf = novel05SnpsDf[order(-novel05SnpsDf$all_af1),]
    novel05SnpsNotrfDf = novel05SnpsDf[is.na(novel05SnpsDf$simple_repeat),]

    nrow(novel05SnpsDf)
    nrow(novel05SnpsNotrfDf)

    colnames(novel05SnpsDf)
    #head(novel0501SnpsDf$repeat_masker)
    #head(novel0501SnpsDf$simple_repeat)
    #unique(novel0501SnpsDf$simple_repeat)
    #head(novel0501SnpsDf$nested_repeat)

    # Regional breakdown of novel SNPs (AF >= 0.05)
    novel05SnpsRegGrpDf = sqldf("SELECT region, COUNT(*) AS cnt FROM novel05SnpsDf GROUP BY region")

    # Mutation spectra of novel SNPs (AF >= 0.05)
    novel05SnpsDf2 = novel05SnpsDf
    novel05SnpsDf2$new_ref = novel05SnpsDf2$ref
    novel05SnpsDf2$new_alt = novel05SnpsDf2$alt

    novel05SnpsDf2[novel05SnpsDf2$ref == "G", "new_ref"] = "C"
    novel05SnpsDf2[novel05SnpsDf2$ref == "G" & novel05SnpsDf2$alt == "A", "new_alt"] = "T"
    novel05SnpsDf2[novel05SnpsDf2$ref == "G" & novel05SnpsDf2$alt == "C", "new_alt"] = "G"
    novel05SnpsDf2[novel05SnpsDf2$ref == "G" & novel05SnpsDf2$alt == "T", "new_alt"] = "A"

    novel05SnpsDf2[novel05SnpsDf2$ref == "A", "new_ref"] = "T"
    novel05SnpsDf2[novel05SnpsDf2$ref == "A" & novel05SnpsDf2$alt == "G", "new_alt"] = "C"
    novel05SnpsDf2[novel05SnpsDf2$ref == "A" & novel05SnpsDf2$alt == "C", "new_alt"] = "G"
    novel05SnpsDf2[novel05SnpsDf2$ref == "A" & novel05SnpsDf2$alt == "T", "new_alt"] = "A"

    novel05SnpsDf2$Substitution_Class = sprintf("%s>%s", novel05SnpsDf2$new_ref, novel05SnpsDf2$new_alt)
    novel05SnpsDf2$Substitution_Class = factor(novel05SnpsDf2$Substitution_Class, level = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))

    novel05SnpsSubGrpDf = sqldf('SELECT Substitution_Class, COUNT(*) AS Count FROM novel05SnpsDf2 GROUP BY Substitution_Class')
    novel05SnpsSubGrpDf$Type = "Novel SNPs"

    # PhastCons of novel SNPs
    novel05SnpsDf2$phast_cons_elements46way_score = as.numeric(gsub(".*lod=(\\d+).*", "\\1", novel05SnpsDf2$phast_cons_elements46way))
    novel05SnpsDf2$phast_cons_elements46way_placental_score = as.numeric(gsub(".*lod=(\\d+).*", "\\1", novel05SnpsDf2$phast_cons_elements46way_placental))
    novel05SnpsDf2$phast_cons_elements46way_primates_score = as.numeric(gsub(".*lod=(\\d+).*", "\\1", novel05SnpsDf2$phast_cons_elements46way_primates))
    boxplot(novel05SnpsDf2$phast_cons_elements46way_score)
    head(novel05SnpsDf2$phast_cons_elements46way_score)
    summary(novel05SnpsDf2$phast_cons_elements46way_score)

    # Plot boxplots for PhastCons distributions
    phastConsDf = rbind(data.frame(type = "Known SNPs", knownSnpsDf2[,c("chrom", "pos", "phast_cons_elements46way_score")]),
                        data.frame(type = "Over-represented known SNPs", orSnpsDf2[,c("chrom", "pos", "phast_cons_elements46way_score")]),
                        data.frame(type = "Novel SNPs", novel05SnpsDf2[,c("chrom", "pos", "phast_cons_elements46way_score")]))

    phastConsPlacentalDf = rbind(data.frame(type = "Known SNPs", knownSnpsDf2[,c("chrom", "pos", "phast_cons_elements46way_placental_score")]),
                                 data.frame(type = "Over-represented known SNPs", orSnpsDf2[,c("chrom", "pos", "phast_cons_elements46way_placental_score")]),
                                 data.frame(type = "Novel SNPs", novel05SnpsDf2[,c("chrom", "pos", "phast_cons_elements46way_placental_score")]))

    phastConsPrimatesDf = rbind(data.frame(type = "Known SNPs", knownSnpsDf2[,c("chrom", "pos", "phast_cons_elements46way_primates_score")]),
                                data.frame(type = "Over-represented known SNPs", orSnpsDf2[,c("chrom", "pos", "phast_cons_elements46way_primates_score")]),
                                data.frame(type = "Novel SNPs", novel05SnpsDf2[,c("chrom", "pos", "phast_cons_elements46way_primates_score")]))

    #boxplot(phast_cons_elements46way_score~type, data=phastConsDf, ylab="PhastCons 46way", xlab="")
    #boxplot(phast_cons_elements46way_placental_score~type, data=phastConsPlacentalDf, ylab="PhastCons score (46way placental)", xlab="")
    boxplot(phast_cons_elements46way_primates_score~type, data=phastConsPrimatesDf, ylab="PhastCons score (46way primates)", xlab="")

    # Statistical test to compare phastCons scores
    wilcox.test(knownSnpsDf2$phast_cons_elements46way_score, orSnpsDf2$phast_cons_elements46way_score)
    wilcox.test(knownSnpsDf2$phast_cons_elements46way_score, novel05SnpsDf2$phast_cons_elements46way_score)
    wilcox.test(orSnpsDf2$phast_cons_elements46way_score, novel05SnpsDf2$phast_cons_elements46way_score)
    t.test(knownSnpsDf2$phast_cons_elements46way_score, orSnpsDf2$phast_cons_elements46way_score)
    t.test(knownSnpsDf2$phast_cons_elements46way_score, novel05SnpsDf2$phast_cons_elements46way_score)
    t.test(orSnpsDf2$phast_cons_elements46way_score, novel05SnpsDf2$phast_cons_elements46way_score)

    wilcox.test(knownSnpsDf2$phast_cons_elements46way_primates_score, orSnpsDf2$phast_cons_elements46way_primates_score)
    wilcox.test(knownSnpsDf2$phast_cons_elements46way_primates_score, novel05SnpsDf2$phast_cons_elements46way_primates_score)
    wilcox.test(orSnpsDf2$phast_cons_elements46way_primates_score, novel05SnpsDf2$phast_cons_elements46way_primates_score)
    t.test(knownSnpsDf2$phast_cons_elements46way_primates_score, orSnpsDf2$phast_cons_elements46way_primates_score)
    t.test(knownSnpsDf2$phast_cons_elements46way_primates_score, novel05SnpsDf2$phast_cons_elements46way_primates_score)
    t.test(orSnpsDf2$phast_cons_elements46way_primates_score, novel05SnpsDf2$phast_cons_elements46way_primates_score)


    baseFontSize = 20
    breaks = 10**(1:4)
    p = ggplot(phastConsDf, aes(x=type, y=phast_cons_elements46way_score)) +
        geom_boxplot() +
        theme(plot.title   = element_text(size = baseFontSize + 3, face="bold"),
                axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
                axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
                axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 30, hjust=1),
                axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
                strip.text.x = element_text(size = baseFontSize, face="bold"),
                strip.text.y = element_text(size = baseFontSize, face="bold", angle = 0),
                legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
                legend.text  = element_text(size = baseFontSize, family="sans"),
                legend.direction = "horizontal",
                legend.position = "bottom") +
        scale_x_discrete(name="") +
        scale_y_log10("PhastCons score (46way)\n", breaks = breaks, labels = comma(breaks))
        #scale_y_continuous("PhastCons score (46way)\n", labels=comma)
        #guides(fill = guide_legend(title = NULL, title.position = "top", nrow = 2, byrow = TRUE))
    #show(p)
    plotFile = file.path(mutFigDir, "PhastCons_46way.lod.pdf")
    ggsave(plot=p, plotFile, width=6, height=8)

    baseFontSize = 20
    breaks = 10**(1:4)
    p = ggplot(phastConsPrimatesDf, aes(x=type, y=phast_cons_elements46way_primates_score)) +
        geom_boxplot() +
        theme(plot.title   = element_text(size = baseFontSize + 3, face="bold"),
                axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
                axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
                axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 30, hjust=1),
                axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
                strip.text.x = element_text(size = baseFontSize, face="bold"),
                strip.text.y = element_text(size = baseFontSize, face="bold", angle = 0),
                legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
                legend.text  = element_text(size = baseFontSize, family="sans"),
                legend.direction = "horizontal",
                legend.position = "bottom") +
        scale_x_discrete(name="") +
        scale_y_log10("PhastCons score (46way primates)\n", breaks = breaks, labels = comma(breaks))
        #scale_y_continuous("PhastCons score (46way primates)\n", labels=comma)
        #guides(fill = guide_legend(title = NULL, title.position = "top", nrow = 2, byrow = TRUE))
    #show(p)
    plotFile = file.path(mutFigDir, "PhastCons_46way_primates.lod.pdf")
    ggsave(plot=p, plotFile, width=6, height=8)


    # Plot boxplots for CADD distributions
    caddRawDf = rbind(data.frame(type = "Known SNPs", knownSnpsDf2[,c("chrom", "pos", "cadd_raw1")]),
                      data.frame(type = "Over-represented known SNPs", orSnpsDf2[,c("chrom", "pos", "cadd_raw1")]),
                      data.frame(type = "Novel SNPs", novel05SnpsDf2[,c("chrom", "pos", "cadd_raw1")]))

    caddPhredDf = rbind(data.frame(type = "Known SNPs", knownSnpsDf2[,c("chrom", "pos", "cadd_phred1")]),
                        data.frame(type = "Over-represented known SNPs", orSnpsDf2[,c("chrom", "pos", "cadd_phred1")]),
                        data.frame(type = "Novel SNPs", novel05SnpsDf2[,c("chrom", "pos", "cadd_phred1")]))

    boxplot(cadd_raw1~type, data=caddRawDf, ylab="CADD Raw score", xlab="")
    boxplot(cadd_phred1~type, data=caddPhredDf, ylab="CADD Phred score", xlab="")

    # Statistical test to compare phastCons scores
    wilcox.test(knownSnpsDf2$cadd_raw1, orSnpsDf2$cadd_raw1)
    wilcox.test(knownSnpsDf2$cadd_raw1, novel05SnpsDf2$cadd_raw1)
    wilcox.test(orSnpsDf2$cadd_raw1, novel05SnpsDf2$cadd_raw1)
    t.test(knownSnpsDf2$cadd_raw1, orSnpsDf2$cadd_raw1)
    t.test(knownSnpsDf2$cadd_raw1, novel05SnpsDf2$cadd_raw1)
    t.test(orSnpsDf2$cadd_raw1, novel05SnpsDf2$cadd_raw1)

    wilcox.test(knownSnpsDf2$cadd_phred1, orSnpsDf2$cadd_phred1)
    wilcox.test(knownSnpsDf2$cadd_phred1, novel05SnpsDf2$cadd_phred1)
    wilcox.test(orSnpsDf2$cadd_phred1, novel05SnpsDf2$cadd_phred1)
    t.test(knownSnpsDf2$cadd_phred1, orSnpsDf2$cadd_phred1)
    t.test(knownSnpsDf2$cadd_phred1, novel05SnpsDf2$cadd_phred1)
    t.test(orSnpsDf2$cadd_phred1, novel05SnpsDf2$cadd_phred1)

    summary(knownSnpsDf2$cadd_raw1)

    baseFontSize = 20
    p = ggplot(caddRawDf, aes(x=type, y=cadd_raw1)) +
        geom_boxplot() +
        theme(plot.title   = element_text(size = baseFontSize + 3, face="bold"),
                axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
                axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
                axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 30, hjust=1),
                axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
                strip.text.x = element_text(size = baseFontSize, face="bold"),
                strip.text.y = element_text(size = baseFontSize, face="bold", angle = 0),
                legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
                legend.text  = element_text(size = baseFontSize, family="sans"),
                legend.direction = "horizontal",
                legend.position = "bottom") +
        scale_x_discrete(name="") +
        #scale_y_continuous("CADD Raw score\n", labels=comma, limits = c(-10, 10))
        scale_y_continuous("CADD Raw score\n", labels=comma, trans = "log2")
    plotFile = file.path(mutFigDir, "CADD_Raw.2.png")
    ggsave(plot=p, plotFile, width=6, height=8)

    baseFontSize = 20
    p = ggplot(caddPhredDf, aes(x=type, y=cadd_phred1)) +
        geom_violin() +
        theme(plot.title   = element_text(size = baseFontSize + 3, face="bold"),
                axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
                axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
                axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 30, hjust=1),
                axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
                strip.text.x = element_text(size = baseFontSize, face="bold"),
                strip.text.y = element_text(size = baseFontSize, face="bold", angle = 0),
                legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
                legend.text  = element_text(size = baseFontSize, family="sans"),
                legend.direction = "horizontal",
                legend.position = "bottom") +
        scale_x_discrete(name="") +
        scale_y_continuous("CADD Phred-scaled score\n", labels=comma)
    plotFile = file.path(mutFigDir, "CADD_Phred.violin.pdf")
    ggsave(plot=p, plotFile, width=6, height=8)



    # Plot mutation spectra
    mutSptraDf = rbind(knownSnpsSubGrpDf, orSnpsSubGrpDf, novel05SnpsSubGrpDf)
    mutSptraDf = subset(mutSptraDf, !is.na(Substitution_Class))
    mutSptraDf$Type = factor(mutSptraDf$Type, levels = c("Known SNPs", "Over-represented known SNPs", "Novel SNPs"))

    baseFontSize = 18
    p = ggplot(mutSptraDf, aes(x=Type, y=Count, fill = Substitution_Class)) +
    geom_bar(stat="identity", position="fill") +
    theme(plot.title   = element_text(size = baseFontSize + 3, face="bold"),
            axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
            axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
            axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 30, hjust=1),
            axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
            strip.text.x = element_text(size = baseFontSize, face="bold"),
            strip.text.y = element_text(size = baseFontSize, face="bold", angle = 0),
            legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
            legend.text  = element_text(size = baseFontSize, family="sans"),
            legend.direction = "horizontal",
            legend.position = "bottom") +
    scale_x_discrete(name="") +
    scale_y_continuous("Proportion\n", labels=comma) +
    scale_fill_brewer(palette="Set1") +
    guides(fill = guide_legend(title = NULL, title.position = "top", nrow = 2, byrow = TRUE))

    #show(p)
    plotFile = file.path(mutFigDir, "Mutation_Spectra.ratio.pdf")
    ggsave(plot=p, plotFile, width=5, height=8)

    baseFontSize = 20
    p = ggplot(mutSptraDf, aes(x=Type, y=Count, fill = Substitution_Class)) +
    geom_bar(stat="identity") +
    theme(plot.title   = element_text(size = baseFontSize + 3, face="bold"),
            axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
            axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
            axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 30, hjust=1),
            axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
            strip.text.x = element_text(size = baseFontSize, face="bold"),
            strip.text.y = element_text(size = baseFontSize, face="bold", angle = 0),
            legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
            legend.text  = element_text(size = baseFontSize, family="sans"),
            legend.direction = "horizontal",
            legend.position = "bottom") +
    scale_x_discrete(name="") +
    scale_y_continuous("Count\n", labels=comma) +
    scale_fill_brewer(palette="Set1") +
    guides(fill = guide_legend(title = NULL, title.position = "top", nrow = 2, byrow = TRUE))

    #show(p)
    plotFile = file.path(mutFigDir, "Mutation_Spectra.count.pdf")
    ggsave(plot=p, plotFile, width=5, height=8)


    novel05SnpsSampleIds = grep("TCGA", colnames(novel05SnpsNotrfDf), value = T)
    novel05SnpsPatientIds = as.vector(unlist(sapply(novel05SnpsSampleIds, function(x) { paste(strsplit(x, "_", fixed = T)[[1]][1:3], collapse = "_") })))
    
    # Patient IDs shared by MAF, SCNA and germline calls
    #sharedPatientIds = intersect(mafPatientIds, novel05SnpsPatientIds)
    sharedPatientIds = intersect(intersect(mafPatientIds, novel05SnpsPatientIds), scnaPatientIds)

    # Exonic sites
    novel05SnpsExonicDf = subset(novel05SnpsNotrfDf, region == "exonic")
    novel05SnpsExonicNonsynDf = subset(novel05SnpsNotrfDf, region == "exonic" & exonic_vartype != "synonymous SNV")
    novel05SnpsExonicNonsynCdDf = subset(novel05SnpsNotrfDf, region == "exonic" & exonic_vartype != "synonymous SNV" & feature %in% cancerGenes)

    #orSnpsDf = orSnpsDf[order(-orSnpsDf$all_af1),]
    #orSnpsExonicDf = subset(orSnpsDf, region == "exonic")
    #orSnpsExonicNonsynDf = subset(orSnpsDf, region == "exonic" & exonic_vartype != "synonymous SNV")
    #orSnpsExonicNonsynCdDf = subset(orSnpsDf, region == "exonic" & exonic_vartype != "synonymous SNV" & feature %in% cancerGenes)

    nrow(novel05SnpsExonicDf)
    nrow(novel05SnpsExonicNonsynDf)
    nrow(novel05SnpsExonicNonsynCdDf)

    novel05SnpsExonicFile = file.path(tableDir,  sprintf("GCC-Germline-%s-Exonic.txt", canType))
    write.table(novel05SnpsExonicDf, file=novel05SnpsExonicFile, sep="\t", row.names=F, col.names=T, quote=FALSE)

    novel05SnpsExonicNonsynFile = file.path(tableDir,  sprintf("GCC-Germline-%s-Exonic-Nonsyn.txt", canType))
    write.table(novel05SnpsExonicNonsynDf, file=novel05SnpsExonicNonsynFile, sep="\t", row.names=F, col.names=T, quote=FALSE)

    novel05SnpsExonicNonsynCdFile = file.path(tableDir,  sprintf("GCC-Germline-%s-Exonic-Nonsyn-Cd.txt", canType))
    write.table(novel05SnpsExonicNonsynCdDf, file=novel05SnpsExonicNonsynCdFile, sep="\t", row.names=F, col.names=T, quote=FALSE)

    # Splicing sites
    novel05SnpsSpliceDf = subset(novel05SnpsNotrfDf, region == "splicing")
    novel05SnpsSpliceCdDf = subset(novel05SnpsNotrfDf, region == "splicing" & feature %in% cancerGenes)
    nrow(novel05SnpsSpliceDf)
    nrow(novel05SnpsSpliceCdDf)
    head(novel05SnpsSpliceDf)
    head(novel05SnpsSpliceCdDf)
    novel05SnpsSpliceFile = file.path(tableDir,  sprintf("GCC-Germline-%s-Splice.txt", canType))
    write.table(novel05SnpsSpliceDf, file=novel05SnpsSpliceFile, sep="\t", row.names=F, col.names=T, quote=FALSE)

    ## Obviously functional germline variants hitting known cancer driver genes
    #novel05SnpsObvFuncDf = rbind(novel05SnpsExonicNonsynCdDf, novel05SnpsSpliceCdDf)
    #novel05SnpsObvFuncDf = novel05SnpsObvFuncDf[order(-novel05SnpsObvFuncDf$all_af1),]
    #novel05SnpsObvFuncFile = file.path(tableDir,  sprintf("GCC-Germline-%s-Exonic_Splicing_SNPs_in_Cancer_Drivers.txt", canType))
    #write.table(novel05SnpsObvFuncDf, file=novel05SnpsObvFuncFile, sep="\t", row.names=F, col.names=T, quote=FALSE)

    # Non-coding but potentially functional
    novel05SnpsNoncoDf = subset(novel05SnpsNotrfDf, region != "exonic" & region != "splicing" & region != "ncRNA_exonic" & region != "ncRNA_splicing")
    novel05SnpsNoncoDf = novel05SnpsNoncoDf[order(-novel05SnpsNoncoDf$all_af1),]
    nrow(novel05SnpsNoncoDf)
    head(novel05SnpsNoncoDf)
    unique(novel05SnpsNoncoDf$region)
    summary(novel05SnpsNoncoDf$all_af1)

    # Test mRNA expression changes of nearby/associated gene(s)
    novel05SnpsNoncoSigExpDf = data.frame()

    for (i in 1:nrow(novel05SnpsNoncoDf)) {
        print(i)
        r = novel05SnpsNoncoDf[i,]

        mutSampleIds = colnames(r)[grepl("'\\S{1}/\\S{1}'", r, perl = TRUE) & 
                                   !grepl("0/0", r, perl = TRUE) & 
                                   !grepl("\\./\\.", r, perl = TRUE)]
        mutPatientIds = gsub("(TCGA_\\S{2}_\\S{4})_\\S{3}", "\\1", mutSampleIds)
        wdtSampleIds = colnames(r)[grepl("0/0", r, perl = TRUE)]
        wdtPatientIds = gsub("(TCGA_\\S{2}_\\S{4})_\\S{3}", "\\1", wdtSampleIds)

        nearbyGeneNames = strsplit(gsub("\\(.*?\\)", "", r$feature), ",")[[1]]
        nearbyGeneDists = as.numeric(gsub(".*dist=(\\d+).*", "\\1", strsplit(r$feature, ",")[[1]]))
        corrDhsDhsGeneNames = strsplit(gsub("\\(.*?\\)", "", gsub("Name=", "", r$distal_dhs_to_promoter_dhs)), ",")[[1]]
        corrDhsRnaGeneNames = strsplit(gsub("\\(.*?\\)", "", gsub("Name=", "", r$dhs_to_gene_expression)), ",")[[1]]
        geneNames = unique(na.omit(c(nearbyGeneNames, corrDhsDhsGeneNames, corrDhsRnaGeneNames)))

        for (geneName in geneNames) {
            if (geneName == "CAAP1") geneName = "C9orf82"
            varStem = sprintf("%s_%d_%s_%s_%s", r$chrom, r$pos, r$ref, r$alt, geneName)
            figStem = sprintf("chr%s:%d:%s>%s", r$chrom, r$pos, r$ref, r$alt)
            selGeneExps = geneExpMergedDf[grep(paste(geneName, "|", sep=""), rownames(geneExpMergedDf), fixed = T),]
            if (nrow(selGeneExps) != 1) {
                print(sprintf("No expression data for %s", geneName))
                #novel05SnpsNoncoSigExpDf = rbind(novel05SnpsNoncoSigExpDf,
                                                       #data.frame(r, 
                                                                  #gene = geneName, 
                                                                  #mean_mutant_gene_expression = NA,
                                                                  #mean_wildtype_gene_expression = NA,
                                                                  #mean_total_gene_expression = NA,
                                                                  #median_mutant_gene_expression = NA,
                                                                  #median_wildtype_gene_expression = NA,
                                                                  #median_total_gene_expression = NA,
                                                                  #mean_mutant_gene_expression_log2_fold_change = NA,
                                                                  #mean_wildtype_gene_expression_log2_fold_change = NA,
                                                                  #mean_total_gene_expression_log2_fold_change = NA,
                                                                  #median_mutant_gene_expression_log2_fold_change = NA,
                                                                  #median_wildtype_gene_expression_log2_fold_change = NA,
                                                                  #median_total_gene_expression_log2_fold_change = NA,
                                                                  #wt_pvalue_gene_expression_mutant_vs_wildtype = NA,
                                                                  #tt_pvalue_gene_expression_mutant_vs_wildtype = NA,
                                                                  #wt_pvalue_gene_expression_mutant_vs_total = NA,
                                                                  #tt_pvalue_gene_expression_mutant_vs_total = NA,
                                                                  #wt_pvalue_gene_expression_log2_fold_change_mutant_vs_wildtype = NA,
                                                                  #tt_pvalue_gene_expression_log2_fold_change_mutant_vs_wildtype = NA,
                                                                  #wt_pvalue_gene_expression_log2_fold_change_mutant_vs_total = NA,
                                                                  #tt_pvalue_gene_expression_log2_fold_change_mutant_vs_total = NA))
                next
            }
      
            mutGeneExps = selGeneExps[, which(gsub("(TCGA_\\S{2}_\\S{4})_\\S{3}", "\\1", colnames(selGeneExps)) %in% mutPatientIds)]
            wdtGeneExps = selGeneExps[, which(gsub("(TCGA_\\S{2}_\\S{4})_\\S{3}", "\\1", colnames(selGeneExps)) %in% wdtPatientIds)]

            ## RSEM
            selGeneExpMutVsWdtDf = data.frame()
            selGeneExpMutVsWdtDf = rbind(selGeneExpMutVsWdtDf, data.frame(type=rep(figStem, ncol(mutGeneExps)), rsem = as.numeric(mutGeneExps[1,])))
            selGeneExpMutVsWdtDf = rbind(selGeneExpMutVsWdtDf, data.frame(type=rep("WT", ncol(wdtGeneExps)), rsem = as.numeric(wdtGeneExps[1,])))

            wtPvalueGeneExpMutVsWdt = my.wilcox.test.p.value(as.numeric(mutGeneExps[1,]), as.numeric(wdtGeneExps[1,]))
            ttPvalueGeneExpMutVsWdt = my.t.test.p.value(as.numeric(mutGeneExps[1,]), as.numeric(wdtGeneExps[1,]))

            pvalueCondGeneExpMutVsWdt = (wtPvalueGeneExpMutVsWdt < 0.05 | ttPvalueGeneExpMutVsWdt < 0.05)

            if (!is.na(pvalueCondGeneExpMutVsWdt) & pvalueCondGeneExpMutVsWdt == TRUE) {
                mutVsWdtExpPlot = ggplot(selGeneExpMutVsWdtDf, aes(factor(type), rsem)) +
                ggtitle(sprintf("%s\n", geneName)) +
                geom_boxplot() +
                theme(plot.title   = element_text(size = baseFontSize + 2, face="bold"),
                      axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
                      axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
                      axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 0, hjust = NULL),
                      axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(color = "black", fill="white"),
                      legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
                      legend.text  = element_text(size = baseFontSize, family="sans"),
                      legend.direction = "horizontal",
                      legend.position = "bottom") +
                                   scale_x_discrete(name=sprintf("\nWilcoxon signed-rank test: %f\nStudent's t-test: %f", wtPvalueGeneExpMutVsWdt, ttPvalueGeneExpMutVsWdt)) +
                                   scale_y_continuous(name = "RSEM\n")
                                   mutVsWdtExpPlotFile = file.path(expFigDir, sprintf("%s-MutVsWdt-RSEM.pdf", varStem))
                                   ggsave(filename = mutVsWdtExpPlotFile, plot = mutVsWdtExpPlot, width = 4, height = 5)
            }

            #
            selGeneExpMutVsTotDf = data.frame()
            selGeneExpMutVsTotDf = rbind(selGeneExpMutVsTotDf, data.frame(type=rep(figStem, ncol(mutGeneExps)), log2FoldChange = as.numeric(mutGeneExps[1,])))
            selGeneExpMutVsTotDf = rbind(selGeneExpMutVsTotDf, data.frame(type=rep("Total", ncol(selGeneExps)), log2FoldChange = as.numeric(selGeneExps[1,])))

            wtPvalueGeneExpMutVsTot = my.wilcox.test.p.value(as.numeric(mutGeneExps[1,]), as.numeric(selGeneExps[1,]))
            ttPvalueGeneExpMutVsTot = my.t.test.p.value(as.numeric(mutGeneExps[1,]), as.numeric(selGeneExps[1,]))

            pvalueCondGeneExpMutVsTot = (wtPvalueGeneExpMutVsTot < 0.05 | ttPvalueGeneExpMutVsTot < 0.05)

            if (!is.na(pvalueCondGeneExpMutVsTot) & pvalueCondGeneExpMutVsTot == TRUE) {
                mutVsTotExpPlot = ggplot(selGeneExpMutVsTotDf, aes(factor(type), log2FoldChange)) +
                ggtitle(sprintf("%s\n", geneName)) +
                geom_boxplot() +
                theme(plot.title   = element_text(size = baseFontSize + 2, face="bold"),
                      axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
                      axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
                      axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 0, hjust = NULL),
                      axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(color = "black", fill="white"),
                      legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
                      legend.text  = element_text(size = baseFontSize, family="sans"),
                      legend.direction = "horizontal",
                      legend.position = "bottom") +
                                   scale_x_discrete(name=sprintf("\nWilcoxon signed-rank test: %f\nStudent's t-test: %f", wtPvalueGeneExpMutVsTot, ttPvalueGeneExpMutVsTot)) +
                                   scale_y_continuous(name = "RSEM\n")
                                   mutVsTotExpPlotFile = file.path(expFigDir, sprintf("%s-MutVsTot-RSEM.pdf", varStem))
                                   ggsave(filename = mutVsTotExpPlotFile, plot = mutVsTotExpPlot, width = 4, height = 5)
            }

            ## Log2 Fold Change of RSEM
            selGeneExpLog2Fcs = geneExpLog2FcPrimDf[grep(paste(geneName, "|", sep=""), rownames(geneExpLog2FcPrimDf), fixed = T),]
            mutGeneExpLog2Fcs = selGeneExpLog2Fcs[, which(gsub("(TCGA_\\S{2}_\\S{4})_\\S{3}", "\\1", colnames(selGeneExpLog2Fcs)) %in% mutPatientIds)]
            wdtGeneExpLog2Fcs = selGeneExpLog2Fcs[, which(gsub("(TCGA_\\S{2}_\\S{4})_\\S{3}", "\\1", colnames(selGeneExpLog2Fcs)) %in% wdtPatientIds)]

            selGeneExpLog2FcMutVsWdtDf = data.frame()
            selGeneExpLog2FcMutVsWdtDf = rbind(selGeneExpLog2FcMutVsWdtDf, 
                                               data.frame(type=rep(figStem, ncol(mutGeneExpLog2Fcs)), log2FoldChange = as.numeric(mutGeneExpLog2Fcs[1,])))
            selGeneExpLog2FcMutVsWdtDf = rbind(selGeneExpLog2FcMutVsWdtDf, 
                                               data.frame(type=rep("WT", ncol(wdtGeneExpLog2Fcs)), log2FoldChange = as.numeric(wdtGeneExpLog2Fcs[1,])))

            wtPvalueGeneExpLog2FcMutVsWdt = my.wilcox.test.p.value(as.numeric(mutGeneExpLog2Fcs[1,]), as.numeric(wdtGeneExpLog2Fcs[1,]))
            ttPvalueGeneExpLog2FcMutVsWdt = my.t.test.p.value(as.numeric(mutGeneExpLog2Fcs[1,]), as.numeric(wdtGeneExpLog2Fcs[1,]))

            pvalueCondGeneExpLog2FcMutVsWdt = (wtPvalueGeneExpLog2FcMutVsWdt < 0.05 | ttPvalueGeneExpLog2FcMutVsWdt < 0.05)

            if (!is.na(pvalueCondGeneExpLog2FcMutVsWdt) & pvalueCondGeneExpLog2FcMutVsWdt == TRUE) {
                mutVsWdtExpLog2FcPlot = ggplot(selGeneExpLog2FcMutVsWdtDf, aes(factor(type), log2FoldChange)) +
                ggtitle(sprintf("%s\n", geneName)) +
                geom_boxplot() +
                theme(plot.title   = element_text(size = baseFontSize + 2, face="bold"),
                      axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
                      axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
                      axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 0, hjust = NULL),
                      axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(color = "black", fill="white"),
                      legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
                      legend.text  = element_text(size = baseFontSize, family="sans"),
                      legend.direction = "horizontal",
                      legend.position = "bottom") +
                                   scale_x_discrete(name=sprintf("\nWilcoxon signed-rank test: %f\nStudent's t-test: %f", wtPvalueGeneExpLog2FcMutVsWdt, ttPvalueGeneExpLog2FcMutVsWdt)) +
                                   scale_y_continuous(name = "Log2 fold change of RSEM\n")
                                   mutVsWdtExpLog2FcPlotFile = file.path(expFigDir, sprintf("%s-MutVsWdt-Log2FoldChangeRSEM.pdf", varStem))
                                   ggsave(filename = mutVsWdtExpLog2FcPlotFile, plot = mutVsWdtExpLog2FcPlot, width = 4, height = 5)
            }

            #
            selGeneExpLog2FcMutVsTotDf = data.frame()
            selGeneExpLog2FcMutVsTotDf = rbind(selGeneExpLog2FcMutVsTotDf, data.frame(type=rep(figStem, ncol(mutGeneExpLog2Fcs)), log2FoldChange = as.numeric(mutGeneExpLog2Fcs[1,])))
            selGeneExpLog2FcMutVsTotDf = rbind(selGeneExpLog2FcMutVsTotDf, data.frame(type=rep("Total", ncol(selGeneExpLog2Fcs)), log2FoldChange = as.numeric(selGeneExpLog2Fcs[1,])))

            wtPvalueGeneExpLog2FcMutVsTot = my.wilcox.test.p.value(as.numeric(mutGeneExpLog2Fcs[1,]), as.numeric(selGeneExpLog2Fcs[1,]))
            ttPvalueGeneExpLog2FcMutVsTot = my.t.test.p.value(as.numeric(mutGeneExpLog2Fcs[1,]), as.numeric(selGeneExpLog2Fcs[1,]))

            pvalueCondGeneExpLog2FcMutVsTot = (wtPvalueGeneExpLog2FcMutVsTot < 0.05 | ttPvalueGeneExpLog2FcMutVsTot < 0.05)

            if (!is.na(pvalueCondGeneExpLog2FcMutVsTot) & pvalueCondGeneExpLog2FcMutVsTot == TRUE) {
                mutVsTotExpLog2FcPlot = ggplot(selGeneExpLog2FcMutVsTotDf, aes(factor(type), log2FoldChange)) +
                ggtitle(sprintf("%s\n", geneName)) +
                geom_boxplot() +
                theme(plot.title   = element_text(size = baseFontSize + 2, face="bold"),
                      axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
                      axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
                      axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 0, hjust = NULL),
                      axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(color = "black", fill="white"),
                      legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
                      legend.text  = element_text(size = baseFontSize, family="sans"),
                      legend.direction = "horizontal",
                      legend.position = "bottom") +
                                   scale_x_discrete(name=sprintf("\nWilcoxon signed-rank test: %f\nStudent's t-test: %f", 
                                                                 wtPvalueGeneExpLog2FcMutVsTot, 
                                                                 ttPvalueGeneExpLog2FcMutVsTot)) +
                                   scale_y_continuous(name = "Log2 fold change of RSEM\n")
                                   mutVsTotExpLog2FcPlotFile = file.path(expFigDir, sprintf("%s-MutVsTot-Log2FoldChangeRSEM.pdf", varStem))
                                   ggsave(filename = mutVsTotExpLog2FcPlotFile, plot = mutVsTotExpLog2FcPlot, width = 4, height = 5)
            }

            if (!is.na(pvalueCondGeneExpMutVsWdt) & pvalueCondGeneExpMutVsWdt == TRUE) {
                upstreamCgSnpsSigExpDf = rbind(upstreamCgSnpsSigExpDf,
                                                 data.frame(r, 
                                                            gene = geneName, 
                                                            mean_mutant_gene_expression = mean(as.numeric(mutGeneExps)),
                                                            mean_wildtype_gene_expression = mean(as.numeric(wdtGeneExps)),
                                                            mean_total_gene_expression = mean(as.numeric(selGeneExps)),
                                                            median_mutant_gene_expression = median(as.numeric(mutGeneExps)),
                                                            median_wildtype_gene_expression = median(as.numeric(wdtGeneExps)),
                                                            median_total_gene_expression = median(as.numeric(selGeneExps)),
                                                            mean_mutant_gene_expression_log2_fold_change = mean(as.numeric(mutGeneExpLog2Fcs)),
                                                            mean_wildtype_gene_expression_log2_fold_change = mean(as.numeric(wdtGeneExpLog2Fcs)),
                                                            mean_total_gene_expression_log2_fold_change = mean(as.numeric(selGeneExpLog2Fcs)),
                                                            median_mutant_gene_expression_log2_fold_change = median(as.numeric(mutGeneExpLog2Fcs)),
                                                            median_wildtype_gene_expression_log2_fold_change = median(as.numeric(wdtGeneExpLog2Fcs)),
                                                            median_total_gene_expression_log2_fold_change = median(as.numeric(selGeneExpLog2Fcs)),
                                                            wt_pvalue_gene_expression_mutant_vs_wildtype = wtPvalueGeneExpMutVsWdt,
                                                            tt_pvalue_gene_expression_mutant_vs_wildtype = ttPvalueGeneExpMutVsWdt,
                                                            wt_pvalue_gene_expression_mutant_vs_total = wtPvalueGeneExpMutVsTot,
                                                            tt_pvalue_gene_expression_mutant_vs_total = ttPvalueGeneExpMutVsTot,
                                                            wt_pvalue_gene_expression_log2_fold_change_mutant_vs_wildtype = wtPvalueGeneExpLog2FcMutVsWdt,
                                                            tt_pvalue_gene_expression_log2_fold_change_mutant_vs_wildtype = ttPvalueGeneExpLog2FcMutVsWdt,
                                                            wt_pvalue_gene_expression_log2_fold_change_mutant_vs_total = wtPvalueGeneExpLog2FcMutVsTot,
                                                            tt_pvalue_gene_expression_log2_fold_change_mutant_vs_total = ttPvalueGeneExpLog2FcMutVsTot))
            }
        }
    }

    upstreamCgSnpsSigExpDf2 = subset(upstreamCgSnpsSigExpDf,
                                       (!is.na(wt_pvalue_gene_expression_mutant_vs_wildtype) & wt_pvalue_gene_expression_mutant_vs_wildtype <= 0.05) |
                                       (!is.na(tt_pvalue_gene_expression_mutant_vs_wildtype) & tt_pvalue_gene_expression_mutant_vs_wildtype <= 0.05))

    nrow(upstreamCgSnpsSigExpDf2)
    head(upstreamCgSnpsSigExpDf2)



    # Upstream
    novel05SnpsUpstreamDf = subset(novel05SnpsDf, region == "upstream")
    novel05SnpsUpstreamPrex2Df = subset(novel05SnpsUpstreamDf, feature == "PREX2")
    nrow(novel05SnpsUpstreamDf)
    head(novel05SnpsUpstreamDf, 20)
    novel05SnpsUpstreamFile = file.path(tableDir,  sprintf("GCC-Germline-%s-Upstream.txt", canType))
    write.table(novel05SnpsUpstreamDf, file=novel05SnpsUpstreamFile, sep="\t", row.names=F, col.names=T, quote=FALSE)

    # Promoter DHS
    novel05SnpsPromoDf = subset(novel05SnpsSrtDf, !is.na(promoter_dhs_master_known) & region != "upstream")
    nrow(novel05SnpsPromoDf)
    head(novel05SnpsPromoDf)
    novel05SnpsPromoterDHSFile = file.path(tableDir,  sprintf("GCC-Germline-%s-PromoterDHS.txt", canType))
    write.table(novel05SnpsPromoterDHSDf, file=novel05SnpsPromoterDHSFile, sep="\t", row.names=F, col.names=T, quote=FALSE)

    # Distal DHS (DHS-DHS)
    novel05SnpsDhsDhsDf = subset(novel05SnpsDf, !is.na(distal_dhs_to_promoter_dhs))
    nrow(novel05SnpsDhsDhsDf)
    head(novel05SnpsDf,2)

    # Distal DHS (DHS-mRNA)
    novel05SnpsDhsExpDf = subset(novel05SnpsDf, !is.na(dhs_to_gene_expression))

    novel05SnpsWithExpDf = data.frame()

    selColsForNonExonicSnps = c("chrom", "pos", "ref", "alt", 
                                "all_af1", "eur_af1", 
                                "region", "feature", "cadd_phred1",
                                "wg_encode_reg_dnase_clustered_v2",
                                "combined_dhs_peak", "distal_dhs_to_promoter_dhs", "dhs_to_gene_expression",
                                "tfbs_cons_sites",
                                "wg_encode_reg_tfbs_clustered_v3",
                                "wg_encode_broad_hmm_gm12878",
                                "repeat_masker",
                                "gwas_catalog")

    ## Construct binary matrices for MAF and germline calls based on shared patient IDs
    # MAF
    mafNonsilDf$patientId = sapply(mafNonsilDf$Tumor_Sample_Barcode,
                                      function(x) { paste(strsplit(x, "-", fixed = T)[[1]][1:3], collapse = "_") })
    mafNonsilGrpDf = sqldf('SELECT Hugo_Symbol AS gene, patientId, COUNT(*) AS mutCnt 
                           FROM mafNonsilDf GROUP BY gene, patientId')
    mafNonsilGrpCastDf = dcast(mafNonsilGrpDf, "gene ~ patientId", value.var = "mutCnt")
    rownames(mafNonsilGrpCastDf) = mafNonsilGrpCastDf$gene
    mafNonsilGrpCastDf = mafNonsilGrpCastDf[,-1]
    mafNonsilGrpCastDf[is.na(mafNonsilGrpCastDf)] = 0
    mafNonsilGrpCastDf[mafNonsilGrpCastDf > 0] = 1 
    mafNonsilGrpCastSharedDf = mafNonsilGrpCastDf[,sharedPatientIds]
    # SCNA
    scnaBinDf = scnaDf[,sharedPatientIds]
    scnaBinDf[scnaBinDf > 0.8 | scnaBinDf < -0.8] = 1
    scnaBinDf[scnaBinDf != 1] = 0
    # Germline SNP
    novel05SnpsGtsDf = novel05SnpsDf[, grep("TCGA", colnames(novel05SnpsDf))]
    rownames(novel05SnpsGtsDf) = gsub(" ", "",
                                       apply(novel05SnpsDf[, c("chrom", "pos", "ref", "alt")], 1, paste, collapse = "_"))
    colnames(novel05SnpsGtsDf) =
    as.vector(unlist(sapply(colnames(novel05SnpsGtsDf), function(x) {
                            paste(strsplit(x, "_", fixed = T)[[1]][1:3],
                                  collapse = "_") })))
    novel05SnpsGtsSharedDf = novel05SnpsGtsDf[, sharedPatientIds]
    novel05SnpsGtsSharedDf[novel05SnpsGtsSharedDf == "'0/0'"] = 0
    novel05SnpsGtsSharedDf[novel05SnpsGtsSharedDf == "'./.'"] = 0
    novel05SnpsGtsSharedDf[novel05SnpsGtsSharedDf != 0] = 1
    novel05SnpsGtsSharedMat = data.matrix(novel05SnpsGtsSharedDf)
    rownames(novel05SnpsGtsSharedMat) = rownames(novel05SnpsGtsSharedDf)
    novel05SnpsGtsSharedDf = as.data.frame(novel05SnpsGtsSharedMat)


    fimoOutDf = data.frame()

    for (j in rownames(novel05SnpsDhsExpDf)) {
        #r = novel05SnpsDf[j,]
        #r = subset(novel05SnpsDf, pos == 54963211) # [AURKA]
        #r = subset(novel05SnpsDf, pos == 68864149) # [PREX2]
        r = subset(novel05SnpsDf, pos == 41740443) # [INHBA]
        #r = subset(novel05SnpsDf, pos == 33108116) # [RBBP4]
        #r = subset(novel05SnpsDf, pos == 1004316) # [TRIP13]
        #r = subset(novel05SnpsDf, pos == 115438806) # [CASP7]
        #r = subset(novel05SnpsDf, pos == 52551362) # [NR4A1]
        #r = subset(novel05SnpsDf, pos == 23822373) # [SLC22A17]
        #r = subset(novel05SnpsDf, pos == 33602540) # [RYR3]
        #r = subset(novel05SnpsDf, pos == 46806235) # [HOXB13]
        #r = subset(novel05SnpsDf, pos == 51255978) # [KLK1]

        if (is.na(r$feature)) {
            print(sprintf("No genes nearby. Skipped"))
            next
        }
        mutSampleIds = colnames(r)[grepl("'\\S{1}/\\S{1}'", r, perl = TRUE) & 
                                   !grepl("0/0", r, perl = TRUE) & 
                                   !grepl("\\./\\.", r, perl = TRUE)]
        mutPatientIds = gsub("(TCGA_\\S{2}_\\S{4})_\\S{3}", "\\1", mutSampleIds)
        wdtSampleIds = colnames(r)[grepl("0/0", r, perl = TRUE)]
        wdtPatientIds = gsub("(TCGA_\\S{2}_\\S{4})_\\S{3}", "\\1", wdtSampleIds)

        ### Mutual exclusivity analysis
        geneNames = strsplit(gsub("\\(.*?\\)", "", r$feature), ",")[[1]]

        calculateJaccardDistance = function(vec1, vec2) {
            c1 = 0; c2 = 0; c3 = 0; c4 = 0
            for (p in 1:length(vec1)) {
                if (vec1[p] == 0 & vec2[p] == 0) {
                    c1 = c1 + 1
                } else if (vec1[p] == 0 & vec2[p] == 1) {
                    c2 = c2 + 1
                } else if (vec1[p] == 1 & vec2[p] == 0) {
                    c3 = c3 + 1
                } else if (vec1[p] == 1 & vec2[p] == 1) {
                    c4 = c4 + 1
                }
            }
            jd = 1 - (c4 / (c2 + c3 + c4))
            return(jd)
        }

        testMutualExclusivityByFisher = function(vec1, vec2) {
            c1 = 0; c2 = 0; c3 = 0; c4 = 0
            for (p in 1:length(vec1)) {
                if (vec1[p] == 0 & vec2[p] == 0) {
                    c1 = c1 + 1
                } else if (vec1[p] == 0 & vec2[p] == 1) {
                    c2 = c2 + 1
                } else if (vec1[p] == 1 & vec2[p] == 0) {
                    c3 = c3 + 1
                } else if (vec1[p] == 1 & vec2[p] == 1) {
                    c4 = c4 + 1
                }
            }
            f = fisher.test(as.table(matrix(c(c1,c2,c3,c4), ncol = 2, byrow = T)))
            print(as.table(matrix(c(c1,c2,c3,c4), ncol = 2, byrow = T)))
            return(f$p.value)
        }

        testMutualExclusivity = function(vec1, vec2, n = 1000) {
            oriJD = calculateJaccardDistance(vec1, vec2)
            nJDs = c()
            for (q in 1:n) {
                nvec1 = sample(vec1)
                nvec2 = sample(vec2)
                nJD = calculateJaccardDistance(nvec1, nvec2)
                nJDs = c(nJDs, nJD)
            }
            pval = length(nJDs[nJDs > oriJD]) / n
            return(pval)
            #f = fisher.test(as.table(matrix(c(c1,c2,c3,c4), ncol = 2, byrow = T)))
            #print(as.table(matrix(c(c1,c2,c3,c4), ncol = 2, byrow = T)))
            #return(f$p.value)
        }

        for (geneName in geneNames) {
            somaticGeneMutsDf = mafNonsilGrpCastSharedDf[rownames(mafNonsilGrpCastSharedDf) == geneName,]
            somaticGeneScnaDf = scnaBinDf[rownames(scnaBinDf) == geneName,]
            somaticGeneVarsDf = somaticGeneMutsDf + somaticGeneScnaDf
            germlineMutsDf = novel05SnpsGtsSharedDf[rownames(novel05SnpsGtsSharedDf) == sprintf("%s_%d_%s_%s", r$chrom, r$pos, r$ref, r$alt),]
            mergedGeneVarsDf = rbind(germlineMutsDf, somaticGeneMutsDf, somaticGeneScnaDf)

            options(scipen=999)
            mergedGeneVarsDf2 = data.frame()
            for (i in 1:nrow(mergedGeneVarsDf)) {
                newRow = 2^(nrow(mergedGeneVarsDf) - i) * mergedGeneVarsDf[i,]
                mergedGeneVarsDf2 = rbind(mergedGeneVarsDf2, newRow)
            }
            sumCols = colSums(mergedGeneVarsDf2)
            sumCols = rev(sort(sumCols))
            mergedGeneVarsDf3 = mergedGeneVarsDf2[,names(sumCols)]
            smgExcFile = file.path(meTableDir,  sprintf("%s_%d_%s_%s-%s.me.txt", r$chrom, r$pos, r$ref, r$alt, geneName))
            write.table(mergedGeneVarsDf3, file=smgExcFile, sep="\t", row.names=T, col.names=T, quote=FALSE)

            testMutualExclusivity(as.vector(somaticGeneMutsDf), as.vector(germlineMutsDf))
            testMutualExclusivity(as.vector(somaticGeneScnaDf), as.vector(germlineMutsDf))
            testMutualExclusivity(as.vector(somaticGeneMutsDf), as.vector(somaticGeneScnaDf))
            testMutualExclusivity(as.vector(somaticGeneVarsDf), as.vector(germlineMutsDf))

            testMutualExclusivityByFisher(as.vector(somaticGeneMutsDf), as.vector(germlineMutsDf))
            testMutualExclusivityByFisher(as.vector(somaticGeneScnaDf), as.vector(germlineMutsDf))
            testMutualExclusivityByFisher(as.vector(somaticGeneMutsDf), as.vector(somaticGeneScnaDf))
            testMutualExclusivityByFisher(as.vector(somaticGeneVarsDf), as.vector(germlineMutsDf))
        }


        ### Survival analysis
        subPatientIds = c(mutPatientIds, wdtPatientIds)
        subClinDf = subset(clinDf, patientId %in% subPatientIds)
        subClinDf$group = ifelse(subClinDf$patientId %in% mutPatientIds, "Mutant", "Wildtype")
        subClinDf$group = factor(subClinDf$group, levels=c("Mutant", "Wildtype"))
        nrow(subClinDf[subClinDf$group == "Mutant",])
        nrow(subClinDf[subClinDf$group == "Wildtype",])
        head(subClinDf)

        #boxplot(daystodeath ~ group, data =  subClinDf, main = "Days to death")
        #boxplot(daystolastfollowup  ~ group, data =  subClinDf, main = "Days to death")

        #survDiff = survdiff(Surv(daystodeath, vitalstatus) ~ group, data = subClinDf)
        #survTest = surv_test(Surv(daystodeath, vitalstatus) ~ group , data = subClinDf, distribution = approximate(B = 10000))
        survDiff = survdiff(Surv(daystolastfollowup, vitalstatus) ~ group, data = subClinDf)
        survTest = surv_test(Surv(daystolastfollowup, vitalstatus) ~ group , data = subClinDf, distribution = approximate(B = 10000))

        survVarStem = sprintf("%s_%d_%s_%s", r$chrom, r$pos, r$ref, r$alt)
        survFigStem = sprintf("chr%s:%d:%s>%s", r$chrom, r$pos, r$ref, r$alt)
        t(r[, selColsForNonExonicSnps])
        length(mutPatientIds)

        survFigFile = file.path(survFigDir, sprintf("%s-MutVsWdt-Survival.pdf", survVarStem))

        pdf(survFigFile)
        plot(survfit(Surv(daystodeath, vitalstatus) ~ group, data = subClinDf, conf.int=0.95, conf.type="log"), 
             lty = c(1,1), col=c("red","black"),
             mark.time = FALSE, ylab = "Probability (survival function)", xlab = "Survival Time in Days",
             main="Kaplan-Meier estimate with 95% confidence bounds")

        legend("topright", legend = c(paste(c("Mutant", "Wildtype"), ": ", survDiff$n[c(1,2)]," cases" ,sep=""),
                                      paste("p-value (log-rank test)", sprintf("%.4f", pchisq(survDiff$chisq,df=1, lower.tail=F)),sep=": ")),
               lty = c(1,1,0), col=c("red","black", "white"), title = survFigStem, bty = "n")

        dev.off() 

        ### Motif analysis
        delta = 50
        seqStart = r$pos - delta
        seqStart = ifelse(seqStart < 0, 0, seqStart)
        seqEnd = r$pos + delta
        stemId = sprintf("%s_%s_%d_%d_%d", canType, r$chrom, r$pos, seqStart, seqEnd)
        fastaName = sprintf("%s.fasta", stemId)
        fastaFile = file.path(memeFastaDir, fastaName)
        seqWt = getSeq(BSgenome.Hsapiens.UCSC.hg19, names=paste("chr", r$chrom, sep=""), start=seqStart, end=seqEnd)
        seqMt = seqWt

        varPos = r$pos - seqStart + 1
        if (r$ref != "-" & r$alt != "-" & nchar(r$ref) == nchar(r$alt)) { #substitution
            varOffset = 0
            varRelEnd = varPos + nchar(r$ref) - 1
            varRef = subseq(seqWt, varPos, varRelEnd)
            if (varRef != DNAString(r$ref)) {
                stop("Wrong offset for reference sequence!", call. = FALSE)
            }
            subseq(seqMt, varPos, varRelEnd) = DNAString(r$alt)
        } else if (r$ref == "-") { #insertion
            varOffset = nchar(r$alt)
            subseq(seqMt, varPos, width = 0) = DNAString(r$alt)
        } else if (r$alt == "-") { #deletion
            varOffset = -nchar(r$ref)
            varRelEnd = varPos + nchar(r$ref) - 1
            subseq(seqMt, varPos, varRelEnd) = DNAString("")
        }

        seqs = c(as.character(seqWt), as.character(seqMt))
        seqNames = c(sprintf("%s_WT", stemId), sprintf("%s_MT", stemId))
        names(seqs) = seqNames
        seqsDs = DNAStringSet(seqs)
        writeXStringSet(seqsDs, filepath = fastaFile, format = "fasta", width = 200)

        motifDir  = file.path(memeFastaDir, stemId)
        dir.create(motifDir, recursive = TRUE, showWarnings = FALSE)

        fimoCmd = sprintf("fimo -oc \"%s\" %s \"%s\"", motifDir, memeMotifs, fastaFile)
        system(fimoCmd, intern = TRUE, ignore.stderr = TRUE)
        fimoOutTxt = file.path(motifDir, "fimo.txt")
        tmpFimoOutDf = read.delim(fimoOutTxt, header = TRUE, as.is = TRUE,
                                  col.names = c("motifJasparId", "seqName", "start", "end", "strand",
                                                "score", "pValue", "qValue", "matchedSeq"))
        if (nrow(tmpFimoOutDf) < 1) {
            next
        }
        tmpFimoOutDf$space = r$chrom
        tmpFimoOutDf$varStart = varPos
        tmpFimoOutDf$varEnd = varPos
        tmpFimoOutDf$seqStart = seqStart
        tmpFimoOutDf$seqEnd = seqEnd
        tmpFimoOutDf$sampleType = sapply(tmpFimoOutDf$seqName, function(x) tail(unlist(strsplit(x, "_")), 1))

        for (q in 1:nrow(tmpFimoOutDf)) {
            o = tmpFimoOutDf[q,]
            tmpFimoOutDf[q, "offsettedStart"] = o$start
            tmpFimoOutDf[q, "offsettedEnd"] = o$end
            if (o$sampleType == "MT" & o$start >= varPos) {
                tmpFimoOutDf[q, "offsettedStart"] = o$start - varOffset
            }
            if (o$sampleType == "MT" & o$end >= varPos) {
                tmpFimoOutDf[q, "offsettedEnd"] = o$end - varOffset
            }
        }

        tmpFimoOutDf = tmpFimoOutDf[order(tmpFimoOutDf$seqName, 
                                          tmpFimoOutDf$offsettedStart, 
                                          tmpFimoOutDf$offsettedEnd),]

        for (q in 1:nrow(tmpFimoOutDf)) {
            fimoHit                         = tmpFimoOutDf[q,]
            tmpFimoOutDf[q, "motifName"]    = memeMotifIdToNamesDf[memeMotifIdToNamesDf$motifJasparId == fimoHit$motifJasparId,]$motifName
            matchedFimoTissueType           = ifelse(fimoHit$sampleType == "WT", "MT", "WT")
            matchedFimoHit                  = subset(tmpFimoOutDf,
                                                     (motifJasparId == fimoHit$motifJasparId &
                                                      offsettedStart == fimoHit$offsettedStart &
                                                      offsettedEnd == fimoHit$offsettedEnd &
                                                      strand == fimoHit$strand &
                                                      sampleType == matchedFimoTissueType))
            if (nrow(matchedFimoHit) == 1) {
                if (fimoHit$sampleType == "WT") {
                    scoreDiff = matchedFimoHit$score - fimoHit$score
                    tmpFimoOutDf[q, "statusChange"] = ifelse(scoreDiff == 0, "unchanged", ifelse(scoreDiff > 0, "strengthened", "weakened"))
                    scoreDiff = 0
                } else {
                    scoreDiff = fimoHit$score - matchedFimoHit$score
                    tmpFimoOutDf[q, "statusChange"] = ifelse(scoreDiff == 0, "unchanged", ifelse(scoreDiff > 0, "strengthened", "weakened"))
                }
                tmpFimoOutDf[q, "scoreDiff"] = scoreDiff
            } else if (nrow(matchedFimoHit) == 0) {
                if (fimoHit$sampleType == "WT") {
                    scoreDiff = -fimoHit$score
                    tmpFimoOutDf[q, "statusChange"] = "destroyed"
                } else {
                    scoreDiff = fimoHit$score
                    tmpFimoOutDf[q, "statusChange"] = "created"
                }
                tmpFimoOutDf[q, "scoreDiff"] = scoreDiff
            } else {
                logerror("%d", q)
                stop("Impossible to have multiple hits of the same TFBS on the same region!")
            }
        }

        # Draw motif assignment plots
        levelHeight = 1

        ## Assign colors for motifs
        nrMotifs = unique(tmpFimoOutDf$motifName)
        if (length(nrMotifs) > 8) {
            colorFunc = colorRampPalette(brewer.pal(8, "Set2"))
            motifColors = colorFunc(length(nrMotifs))
        } else if (length(nrMotifs) > 2) {
            motifColors = rev(brewer.pal(length(nrMotifs), "Set2"))
        } else {
            motifColors = brewer.pal(3, "Set2")
            motifColors = motifColors[1:length(nrMotifs)]
        }
        names(motifColors) = nrMotifs

        ## Draw wildtype mutant profile
        tmpFimoOutWtDf = subset(tmpFimoOutDf, sampleType == "WT")
        tmpFimoOutPlusWtDf = subset(tmpFimoOutWtDf, strand == "+")
        tmpFimoOutPlusWtGr = as(as(tmpFimoOutPlusWtDf, "RangedData"), "GRanges")
        if (length(tmpFimoOutPlusWtGr) > 0) {
            tmpFimoOutPlusWtGr$level = 0
            for (o in 1:length(tmpFimoOutPlusWtGr)) {
                lvl = 1
                repeat {
                    tmpFimoOutPlusWtLvlGr = subset(tmpFimoOutPlusWtGr, level == lvl)
                    if (length(tmpFimoOutPlusWtLvlGr) == 0) {
                        tmpFimoOutPlusWtGr[o,]$level = lvl
                        break
                    } else {
                        if (countOverlaps(tmpFimoOutPlusWtGr[o,], tmpFimoOutPlusWtLvlGr) == 0) {
                            tmpFimoOutPlusWtGr[o,]$level = lvl
                            break
                        } else {
                            lvl = lvl + 1
                        }
                    }
                }
            }
        }
        tmpFimoOutMinusWtDf = subset(tmpFimoOutWtDf, strand == "-")
        tmpFimoOutMinusWtGr = as(as(tmpFimoOutMinusWtDf, "RangedData"), "GRanges")
        if (length(tmpFimoOutMinusWtGr) > 0) {
            tmpFimoOutMinusWtGr$level = 0
            for (o in 1:length(tmpFimoOutMinusWtGr)) {
                lvl = 1
                repeat {
                    tmpFimoOutMinusWtLvlGr = subset(tmpFimoOutMinusWtGr, level == lvl)
                    if (length(tmpFimoOutMinusWtLvlGr) == 0) {
                        tmpFimoOutMinusWtGr[o,]$level = lvl
                        break
                    } else {
                        if (countOverlaps(tmpFimoOutMinusWtGr[o,], tmpFimoOutMinusWtLvlGr) == 0) {
                            tmpFimoOutMinusWtGr[o,]$level = lvl
                            break
                        } else {
                            lvl = lvl + 1
                        }
                    }
                }
            }
        }
        tmpFimoOutWtGr = c(tmpFimoOutPlusWtGr, tmpFimoOutMinusWtGr)

        # Draw mutant motif profile
        tmpFimoOutMtDf = subset(tmpFimoOutDf, sampleType == "MT")
        tmpFimoOutPlusMtDf = subset(tmpFimoOutMtDf, strand == "+")
        tmpFimoOutPlusMtGr = as(as(tmpFimoOutPlusMtDf, "RangedData"), "GRanges")
        if (length(tmpFimoOutPlusMtGr) > 0) {
            tmpFimoOutPlusMtGr$level = 0
            for (o in 1:length(tmpFimoOutPlusMtGr)) {
                lvl = 1
                repeat {
                    tmpFimoOutPlusMtLvlGr = subset(tmpFimoOutPlusMtGr, level == lvl)
                    if (length(tmpFimoOutPlusMtLvlGr) == 0) {
                        tmpFimoOutPlusMtGr[o,]$level = lvl
                        break
                    } else {
                        if (countOverlaps(tmpFimoOutPlusMtGr[o,], tmpFimoOutPlusMtLvlGr) == 0) {
                            tmpFimoOutPlusMtGr[o,]$level = lvl
                            break
                        } else {
                            lvl = lvl + 1
                        }
                    }
                }
            }
        }

        tmpFimoOutMinusMtDf = subset(tmpFimoOutMtDf, strand == "-")
        tmpFimoOutMinusMtGr = as(as(tmpFimoOutMinusMtDf, "RangedData"), "GRanges")
        if (length(tmpFimoOutMinusMtGr) > 0) {
            tmpFimoOutMinusMtGr$level = 0
            for (o in 1:length(tmpFimoOutMinusMtGr)) {
                lvl = 1
                repeat {
                    tmpFimoOutMinusMtLvlGr = subset(tmpFimoOutMinusMtGr, level == lvl)
                    if (length(tmpFimoOutMinusMtLvlGr) == 0) {
                        tmpFimoOutMinusMtGr[o,]$level = lvl
                        break
                    } else {
                        if (countOverlaps(tmpFimoOutMinusMtGr[o,], tmpFimoOutMinusMtLvlGr) == 0) {
                            tmpFimoOutMinusMtGr[o,]$level = lvl
                            break
                        } else {
                            lvl = lvl + 1
                        }
                    }
                }
            }
        }
        tmpFimoOutMtGr = c(tmpFimoOutPlusMtGr, tmpFimoOutMinusMtGr)

        yEnd1 = if(length(tmpFimoOutPlusWtGr) == 0) { yEnd = 1 } else { yEnd = max(tmpFimoOutPlusWtGr$level) + 2 }
        yStart1 = if(length(tmpFimoOutMinusWtGr) == 0) { yStart = 0 } else { yStart = -max(tmpFimoOutMinusWtGr$level) - 2 }
        yEnd2 = if(length(tmpFimoOutPlusMtGr) == 0) { yEnd = 1 } else { yEnd = max(tmpFimoOutPlusMtGr$level) + 2 }
        yStart2 = if(length(tmpFimoOutMinusMtGr) == 0) { yStart = 0 } else { yStart = -max(tmpFimoOutMinusMtGr$level) - 2 }
        yEnd = max(yEnd1, yEnd2)
        yStart = min(yStart1, yStart2)

        #motifWtFigFile = file.path(memeFastaDir, sprintf("%s.zhao2011.pdf", stemId))
        #motifWtFigFile = file.path(memeFastaDir, sprintf("%s.jolma2013.pdf", stemId))
        motifWtFigFile = file.path(memeFastaDir, sprintf("%s.jaspar_core_2014_vertebrates.pdf", stemId))
        pdf(file = motifWtFigFile, width = 10, height = 9)
        if (length(motifColors) < 10) { lheight = 1 } else { lheight = 2 }
        layout(matrix(c(1, 2, 3, 4), 4, 1, byrow = TRUE), heights = c(4.5, 5, 1, lheight))
        par(cex = 1.3, mar = c(0, 4.1, 2, 1))
        plot(c(seqStart, seqEnd), c(yStart, yEnd), type = "n", xlab = NA, ylab = NA, xaxt = 'n', yaxt = 'n')
        mtext("Wildtype", side = 2, font = 2, cex = 1.5, line = 1)
        for (k in 1:length(tmpFimoOutWtGr)) {
            s = ifelse(toString(strand(tmpFimoOutWtGr[k])) == "-", -1, +1)
            e = tmpFimoOutWtGr[k,]
            lt = 0
            if (e$statusChange == "created" | e$statusChange == "destroyed") {
                lt = 1
            } else if (e$statusChange == "strengthened") {
                lt = 2
            } else if (e$statusChange == "weakened") {
                lt = 3
            }
            rect(seqStart + e$offsettedStart, s * (e$level - 0.9),
                 seqStart + e$offsettedEnd, s * e$level,
                 border = "black", lwd = 2,
                 col = motifColors[toString(e$motifName)], lty = lt)
        }
        abline(h = 0, col = "black", lwd = 5)
        abline(v = seqStart + e$varStart, col = "red", lwd = 3, lty = 2)
        par(xpd = NA)
        text(x = seqStart + e$varStart, yEnd + (yEnd - yStart) / 10, sprintf("chr%s:%d:%s>%s", r$chrom, r$pos, r$ref, r$alt), cex = 1)

        par(cex = 1.3, xpd = FALSE, mar = c(4.1, 4.1, 0, 1))
        plot(c(seqStart, seqEnd), c(yStart, yEnd), type = "n", xlab = sprintf("Genomic position of chromosome %s", r$chrom), ylab = NA, yaxt = 'n')
        mtext("Mutant", side = 2, font = 2, cex = 1.5, line = 1)
        for (k in 1:length(tmpFimoOutMtGr)) {
            s = ifelse(toString(strand(tmpFimoOutMtGr[k])) == "-", -1, +1)
            e = tmpFimoOutMtGr[k,]
            lt = 0
            if (e$statusChange == "created" | e$statusChange == "destroyed") {
                lt = 1
            } else if (e$statusChange == "strengthened") {
                lt = 2
            } else if (e$statusChange == "weakened") {
                lt = 3
            }
            rect(seqStart + e$offsettedStart, s * (e$level - 0.9),
                 seqStart + e$offsettedEnd, s * e$level,
                 border = "black", lwd = 2,
                 col = motifColors[toString(e$motifName)], lty = lt)
        }
        abline(h = 0, col = "black", lwd = 5)
        abline(v = seqStart + e$varStart, col = "red", lwd = 3, lty = 2)

        par(cex = 1.3, xpd = FALSE, mar = c(0, 3, 0, .5))
        plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = NA, ylab = NA)
        legend("bottom", title = "Motif status", legend = c("Strengthened", "Weakened", "Created/Destroyed"),
               lwd = 2, lty = c(2, 3, 1), cex = 0.8, col = "black", bty = "n", ncol = 3, pt.cex = 1.5)

        if (length(motifColors) < 5) { ncol = length(motifColors) } else { ncol = 5 }
        par(cex = 1.3, xpd = FALSE, mar = c(0, 3, 0, .5))
        plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = NA, ylab = NA)
        legend("bottom", title = "TF-binding motif", legend = unique(tmpFimoOutDf$motifName),
               pch = 15, cex = 0.8, col = motifColors, bty = "n", ncol = ncol, pt.cex = 1.5)

        dev.off()

        ### Expression analysis
        geneNames = strsplit(gsub("\\(.*?\\)", "", r$feature), ",")[[1]]
        #geneNames = strsplit(gsub("\\(.*?\\)", "", gsub("Name=", "", r$dhs_to_gene_expression)), ",")[[1]]

        for (geneName in geneNames) {
            if (geneName == "CAAP1") geneName = "C9orf82"
            varStem = sprintf("%s_%d_%s_%s_%s", r$chrom, r$pos, r$ref, r$alt, geneName)
            figStem = sprintf("chr%s:%d:%s>%s", r$chrom, r$pos, r$ref, r$alt)
            selGeneExps = geneExpMergedDf[grep(paste(geneName, "|", sep=""), rownames(geneExpMergedDf), fixed = T),]
            if (nrow(selGeneExps) != 1) {
                print(sprintf("No expression data for %s", geneName))
                novel05SnpsWithExpDf = rbind(novel05SnpsWithExpDf,
                                                       data.frame(r, 
                                                                  gene = geneName, 
                                                                  mean_mutant_gene_expression = NA,
                                                                  mean_wildtype_gene_expression = NA,
                                                                  mean_total_gene_expression = NA,
                                                                  median_mutant_gene_expression = NA,
                                                                  median_wildtype_gene_expression = NA,
                                                                  median_total_gene_expression = NA,
                                                                  mean_mutant_gene_expression_log2_fold_change = NA,
                                                                  mean_wildtype_gene_expression_log2_fold_change = NA,
                                                                  mean_total_gene_expression_log2_fold_change = NA,
                                                                  median_mutant_gene_expression_log2_fold_change = NA,
                                                                  median_wildtype_gene_expression_log2_fold_change = NA,
                                                                  median_total_gene_expression_log2_fold_change = NA,
                                                                  wt_pvalue_gene_expression_mutant_vs_wildtype = NA,
                                                                  tt_pvalue_gene_expression_mutant_vs_wildtype = NA,
                                                                  wt_pvalue_gene_expression_mutant_vs_total = NA,
                                                                  tt_pvalue_gene_expression_mutant_vs_total = NA,
                                                                  wt_pvalue_gene_expression_log2_fold_change_mutant_vs_wildtype = NA,
                                                                  tt_pvalue_gene_expression_log2_fold_change_mutant_vs_wildtype = NA,
                                                                  wt_pvalue_gene_expression_log2_fold_change_mutant_vs_total = NA,
                                                                  tt_pvalue_gene_expression_log2_fold_change_mutant_vs_total = NA))
                next
            }
      
            mutGeneExps = selGeneExps[, which(gsub("(TCGA_\\S{2}_\\S{4})_\\S{3}", "\\1", colnames(selGeneExps)) %in% mutPatientIds)]
            wdtGeneExps = selGeneExps[, which(gsub("(TCGA_\\S{2}_\\S{4})_\\S{3}", "\\1", colnames(selGeneExps)) %in% wdtPatientIds)]

            ## RSEM
            selGeneExpMutVsWdtDf = data.frame()
            selGeneExpMutVsWdtDf = rbind(selGeneExpMutVsWdtDf, data.frame(type=rep(figStem, ncol(mutGeneExps)), rsem = as.numeric(mutGeneExps[1,])))
            selGeneExpMutVsWdtDf = rbind(selGeneExpMutVsWdtDf, data.frame(type=rep("WT", ncol(wdtGeneExps)), rsem = as.numeric(wdtGeneExps[1,])))

            wtPvalueGeneExpMutVsWdt = my.wilcox.test.p.value(as.numeric(mutGeneExps[1,]), as.numeric(wdtGeneExps[1,]))
            ttPvalueGeneExpMutVsWdt = my.t.test.p.value(as.numeric(mutGeneExps[1,]), as.numeric(wdtGeneExps[1,]))

            pvalueCondGeneExpMutVsWdt = (wtPvalueGeneExpMutVsWdt < 0.05 | ttPvalueGeneExpMutVsWdt < 0.05)

            if (!is.na(pvalueCondGeneExpMutVsWdt) & pvalueCondGeneExpMutVsWdt == TRUE) {
                mutVsWdtExpPlot = ggplot(selGeneExpMutVsWdtDf, aes(factor(type), rsem)) +
                ggtitle(sprintf("%s\n", geneName)) +
                geom_boxplot() +
                theme(plot.title   = element_text(size = baseFontSize + 2, face="bold"),
                      axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
                      axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
                      axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 0, hjust = NULL),
                      axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(color = "black", fill="white"),
                      legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
                      legend.text  = element_text(size = baseFontSize, family="sans"),
                      legend.direction = "horizontal",
                      legend.position = "bottom") +
                                   scale_x_discrete(name=sprintf("\nWilcoxon signed-rank test: %f\nStudent's t-test: %f", wtPvalueGeneExpMutVsWdt, ttPvalueGeneExpMutVsWdt)) +
                                   scale_y_continuous(name = "RSEM\n")
                                   mutVsWdtExpPlotFile = file.path(expFigDir, sprintf("%s-MutVsWdt-RSEM.pdf", varStem))
                                   ggsave(filename = mutVsWdtExpPlotFile, plot = mutVsWdtExpPlot, width = 4, height = 5)
            }

            #
            selGeneExpMutVsTotDf = data.frame()
            selGeneExpMutVsTotDf = rbind(selGeneExpMutVsTotDf, data.frame(type=rep(figStem, ncol(mutGeneExps)), log2FoldChange = as.numeric(mutGeneExps[1,])))
            selGeneExpMutVsTotDf = rbind(selGeneExpMutVsTotDf, data.frame(type=rep("Total", ncol(selGeneExps)), log2FoldChange = as.numeric(selGeneExps[1,])))

            wtPvalueGeneExpMutVsTot = my.wilcox.test.p.value(as.numeric(mutGeneExps[1,]), as.numeric(selGeneExps[1,]))
            ttPvalueGeneExpMutVsTot = my.t.test.p.value(as.numeric(mutGeneExps[1,]), as.numeric(selGeneExps[1,]))

            pvalueCondGeneExpMutVsTot = (wtPvalueGeneExpMutVsTot < 0.05 | ttPvalueGeneExpMutVsTot < 0.05)

            if (!is.na(pvalueCondGeneExpMutVsTot) & pvalueCondGeneExpMutVsTot == TRUE) {
                mutVsTotExpPlot = ggplot(selGeneExpMutVsTotDf, aes(factor(type), log2FoldChange)) +
                ggtitle(sprintf("%s\n", geneName)) +
                geom_boxplot() +
                theme(plot.title   = element_text(size = baseFontSize + 2, face="bold"),
                      axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
                      axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
                      axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 0, hjust = NULL),
                      axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(color = "black", fill="white"),
                      legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
                      legend.text  = element_text(size = baseFontSize, family="sans"),
                      legend.direction = "horizontal",
                      legend.position = "bottom") +
                                   scale_x_discrete(name=sprintf("\nWilcoxon signed-rank test: %f\nStudent's t-test: %f", wtPvalueGeneExpMutVsTot, ttPvalueGeneExpMutVsTot)) +
                                   scale_y_continuous(name = "RSEM\n")
                                   mutVsTotExpPlotFile = file.path(expFigDir, sprintf("%s-MutVsTot-RSEM.pdf", varStem))
                                   ggsave(filename = mutVsTotExpPlotFile, plot = mutVsTotExpPlot, width = 4, height = 5)
            }

            ## Log2 Fold Change of RSEM
            selGeneExpLog2Fcs = geneExpLog2FcPrimDf[grep(paste(geneName, "|", sep=""), rownames(geneExpLog2FcPrimDf), fixed = T),]
            mutGeneExpLog2Fcs = selGeneExpLog2Fcs[, which(gsub("(TCGA_\\S{2}_\\S{4})_\\S{3}", "\\1", colnames(selGeneExpLog2Fcs)) %in% mutPatientIds)]
            wdtGeneExpLog2Fcs = selGeneExpLog2Fcs[, which(gsub("(TCGA_\\S{2}_\\S{4})_\\S{3}", "\\1", colnames(selGeneExpLog2Fcs)) %in% wdtPatientIds)]

            selGeneExpLog2FcMutVsWdtDf = data.frame()
            selGeneExpLog2FcMutVsWdtDf = rbind(selGeneExpLog2FcMutVsWdtDf, 
                                               data.frame(type=rep(figStem, ncol(mutGeneExpLog2Fcs)), log2FoldChange = as.numeric(mutGeneExpLog2Fcs[1,])))
            selGeneExpLog2FcMutVsWdtDf = rbind(selGeneExpLog2FcMutVsWdtDf, 
                                               data.frame(type=rep("WT", ncol(wdtGeneExpLog2Fcs)), log2FoldChange = as.numeric(wdtGeneExpLog2Fcs[1,])))

            wtPvalueGeneExpLog2FcMutVsWdt = my.wilcox.test.p.value(as.numeric(mutGeneExpLog2Fcs[1,]), as.numeric(wdtGeneExpLog2Fcs[1,]))
            ttPvalueGeneExpLog2FcMutVsWdt = my.t.test.p.value(as.numeric(mutGeneExpLog2Fcs[1,]), as.numeric(wdtGeneExpLog2Fcs[1,]))

            pvalueCondGeneExpLog2FcMutVsWdt = (wtPvalueGeneExpLog2FcMutVsWdt < 0.05 | ttPvalueGeneExpLog2FcMutVsWdt < 0.05)

            if (!is.na(pvalueCondGeneExpLog2FcMutVsWdt) & pvalueCondGeneExpLog2FcMutVsWdt == TRUE) {
                mutVsWdtExpLog2FcPlot = ggplot(selGeneExpLog2FcMutVsWdtDf, aes(factor(type), log2FoldChange)) +
                ggtitle(sprintf("%s\n", geneName)) +
                geom_boxplot() +
                theme(plot.title   = element_text(size = baseFontSize + 2, face="bold"),
                      axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
                      axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
                      axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 0, hjust = NULL),
                      axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(color = "black", fill="white"),
                      legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
                      legend.text  = element_text(size = baseFontSize, family="sans"),
                      legend.direction = "horizontal",
                      legend.position = "bottom") +
                                   scale_x_discrete(name=sprintf("\nWilcoxon signed-rank test: %f\nStudent's t-test: %f", wtPvalueGeneExpLog2FcMutVsWdt, ttPvalueGeneExpLog2FcMutVsWdt)) +
                                   scale_y_continuous(name = "Log2 fold change of RSEM\n")
                                   mutVsWdtExpLog2FcPlotFile = file.path(expFigDir, sprintf("%s-MutVsWdt-Log2FoldChangeRSEM.pdf", varStem))
                                   ggsave(filename = mutVsWdtExpLog2FcPlotFile, plot = mutVsWdtExpLog2FcPlot, width = 4, height = 5)
            }

            #
            selGeneExpLog2FcMutVsTotDf = data.frame()
            selGeneExpLog2FcMutVsTotDf = rbind(selGeneExpLog2FcMutVsTotDf, data.frame(type=rep(figStem, ncol(mutGeneExpLog2Fcs)), log2FoldChange = as.numeric(mutGeneExpLog2Fcs[1,])))
            selGeneExpLog2FcMutVsTotDf = rbind(selGeneExpLog2FcMutVsTotDf, data.frame(type=rep("Total", ncol(selGeneExpLog2Fcs)), log2FoldChange = as.numeric(selGeneExpLog2Fcs[1,])))

            wtPvalueGeneExpLog2FcMutVsTot = my.wilcox.test.p.value(as.numeric(mutGeneExpLog2Fcs[1,]), as.numeric(selGeneExpLog2Fcs[1,]))
            ttPvalueGeneExpLog2FcMutVsTot = my.t.test.p.value(as.numeric(mutGeneExpLog2Fcs[1,]), as.numeric(selGeneExpLog2Fcs[1,]))

            pvalueCondGeneExpLog2FcMutVsTot = (wtPvalueGeneExpLog2FcMutVsTot < 0.05 | ttPvalueGeneExpLog2FcMutVsTot < 0.05)

            if (!is.na(pvalueCondGeneExpLog2FcMutVsTot) & pvalueCondGeneExpLog2FcMutVsTot == TRUE) {
                mutVsTotExpLog2FcPlot = ggplot(selGeneExpLog2FcMutVsTotDf, aes(factor(type), log2FoldChange)) +
                ggtitle(sprintf("%s\n", geneName)) +
                geom_boxplot() +
                theme(plot.title   = element_text(size = baseFontSize + 2, face="bold"),
                      axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
                      axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
                      axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 0, hjust = NULL),
                      axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(color = "black", fill="white"),
                      legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
                      legend.text  = element_text(size = baseFontSize, family="sans"),
                      legend.direction = "horizontal",
                      legend.position = "bottom") +
                                   scale_x_discrete(name=sprintf("\nWilcoxon signed-rank test: %f\nStudent's t-test: %f", 
                                                                 wtPvalueGeneExpLog2FcMutVsTot, 
                                                                 ttPvalueGeneExpLog2FcMutVsTot)) +
                                   scale_y_continuous(name = "Log2 fold change of RSEM\n")
                                   mutVsTotExpLog2FcPlotFile = file.path(expFigDir, sprintf("%s-MutVsTot-Log2FoldChangeRSEM.pdf", varStem))
                                   ggsave(filename = mutVsTotExpLog2FcPlotFile, plot = mutVsTotExpLog2FcPlot, width = 4, height = 5)
            }

            novel05SnpsWithExpDf = rbind(novel05SnpsWithExpDf,
                                                   data.frame(r, 
                                                              gene = geneName, 
                                                              mean_mutant_gene_expression = mean(as.numeric(mutGeneExps)),
                                                              mean_wildtype_gene_expression = mean(as.numeric(wdtGeneExps)),
                                                              mean_total_gene_expression = mean(as.numeric(selGeneExps)),
                                                              median_mutant_gene_expression = median(as.numeric(mutGeneExps)),
                                                              median_wildtype_gene_expression = median(as.numeric(wdtGeneExps)),
                                                              median_total_gene_expression = median(as.numeric(selGeneExps)),
                                                              mean_mutant_gene_expression_log2_fold_change = mean(as.numeric(mutGeneExpLog2Fcs)),
                                                              mean_wildtype_gene_expression_log2_fold_change = mean(as.numeric(wdtGeneExpLog2Fcs)),
                                                              mean_total_gene_expression_log2_fold_change = mean(as.numeric(selGeneExpLog2Fcs)),
                                                              median_mutant_gene_expression_log2_fold_change = median(as.numeric(mutGeneExpLog2Fcs)),
                                                              median_wildtype_gene_expression_log2_fold_change = median(as.numeric(wdtGeneExpLog2Fcs)),
                                                              median_total_gene_expression_log2_fold_change = median(as.numeric(selGeneExpLog2Fcs)),
                                                              wt_pvalue_gene_expression_mutant_vs_wildtype = wtPvalueGeneExpMutVsWdt,
                                                              tt_pvalue_gene_expression_mutant_vs_wildtype = ttPvalueGeneExpMutVsWdt,
                                                              wt_pvalue_gene_expression_mutant_vs_total = wtPvalueGeneExpMutVsTot,
                                                              tt_pvalue_gene_expression_mutant_vs_total = ttPvalueGeneExpMutVsTot,
                                                              wt_pvalue_gene_expression_log2_fold_change_mutant_vs_wildtype = wtPvalueGeneExpLog2FcMutVsWdt,
                                                              tt_pvalue_gene_expression_log2_fold_change_mutant_vs_wildtype = ttPvalueGeneExpLog2FcMutVsWdt,
                                                              wt_pvalue_gene_expression_log2_fold_change_mutant_vs_total = wtPvalueGeneExpLog2FcMutVsTot,
                                                              tt_pvalue_gene_expression_log2_fold_change_mutant_vs_total = ttPvalueGeneExpLog2FcMutVsTot))
        }
    }
    #novel05SnpsWithExpFile = file.path(canVcfDir, "novel05SnpsWithExp.txt")
    #write.table(novel05SnpsWithExpDf, novel05SnpsWithExpFile, sep = "\t", quote = FALSE, row.names = FALSE)
}
