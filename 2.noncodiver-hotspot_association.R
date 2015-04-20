args = commandArgs(TRUE)
annotFile = "/n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/hotspot/chrs/Y/TCGA_16_Cancer_Types.wgs.somatic.chrY.sanitized.scnt.hotspot100.fdr0.05.annotated.00019.txt"
annotFile = args[1]
numCores = as.integer(args[2])

##
## Load libraries
##
require(doMC)
registerDoMC(numCores)

require(sqldf)
require(proxy)
require(gdata)
require(gtools)
require(ggplot2)
require(annotate)
require(snpStats)
require(data.table)
require(Biostrings)
require(org.Hs.eg.db)
require(GenomicRanges)
require(VariantAnnotation)
require(BSgenome.Hsapiens.UCSC.hg19)


##
## Initialize global variables
##
rootDir = "/home/sl279"
baseDir = file.path(rootDir, "BiO/Research/NoncoDiver")
vcfDir = file.path(baseDir, "vcf")
hotspotDir = file.path(vcfDir, "pancan")
hotspotDir = file.path(baseDir, "hotspot")
gdacDir = file.path(baseDir, "gdac")
gdacStdDataDir = file.path(gdacDir, "stddata__2015_02_04")
gdacAnlDataDir = file.path(gdacDir, "analyses__2014_10_17")
figDir = file.path(baseDir, "figure")

chrs = c(1:22, "X", "Y")
chrs = factor(chrs, level=chrs)
cchrs = sapply(chrs, function(x) paste("chr", x, sep=""))
cchrs = factor(cchrs, level=cchrs)
baseFontSize = 15

# Load chrom size information
hg19File = file.path(baseDir, "ucsc/database/hg19.genome")
hg19Df = read.delim(hg19File, header=T, as.is=T)
hg19Df = hg19Df[hg19Df$chrom %in% chrs,]
hg19Df$chrom = factor(hg19Df$chrom, levels=chrs)
hg19Df = hg19Df[order(hg19Df$chrom),]


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

## Load RefSeq ID to gene name mapping table
geneInfoFile = file.path(baseDir, "entrez/gene_info_human.txt")
geneInfoDt = fread(geneInfoFile, header = F)
geneInfoDt = geneInfoDt[, c(2:5), with = F]
setnames(geneInfoDt, c("GeneID", "Symbol", "LocusTag", "Synonyms"))
geneInfoExtDt = geneInfoDt[, list(Synonym = unlist(strsplit(Synonyms, "|", fixed = T), use.names = FALSE)), by = "GeneID,Symbol,LocusTag"]
setkey(geneInfoExtDt, Symbol)


## Load total sample IDs for hotspots
vcfFiles = Sys.glob(file.path(vcfDir, "cancer/*/*.stage4.vcf.gz"))
sampleIdsByCancerLst = list()
sampleIdToCancerDf = data.frame()
for (vcfFile in vcfFiles) {
    print(vcfFile)
    cancerType = basename(dirname(vcfFile))
    sampleIds = samples(scanVcfHeader(vcfFile))
    sampleAbbIds = sapply(sampleIds, function(x) { 
                              elems = strsplit(x, "-", fixed = T)[[1]][1:4]
                              elems[1] = "TCGA"
                              elems[4] = substring(elems[4], 1, 2)
                              paste(elems, collapse = "_") })
    sampleIdsByCancerLst[[cancerType]] = sampleAbbIds
    for (sampleAbbId in sampleAbbIds) {
        sampleIdToCancerDf = rbind(sampleIdToCancerDf, data.frame(sid = sampleAbbId, cancer = cancerType))
    }
}
totSids = as.vector(unlist(sampleIdsByCancerLst))


# Load exome SNV/indel data from GDAC
mafRdata = file.path(gdacAnlDataDir, "mafPanFunc.RData")
load(mafRdata)


## Load mRNA expression data from GDAC
geneExpLog2FoldPanRdata = file.path(gdacStdDataDir, "geneExpLog2FoldPan.RData")
load(geneExpLog2FoldPanRdata)
geneExpLog2FoldPanDt = as.data.table(geneExpLog2FoldPanDf)
setkey(geneExpLog2FoldPanDt, symbol)


## Load SCNA data for filter
cnvByGeneRdata = file.path(gdacStdDataDir, "cnvByGene.RData")
load(cnvByGeneRdata)


## Load hotspot data
annotDf = read.delim(annotFile, header = T, as.is = T)
annotDf$allele = NULL
annotDf$strand = NULL
annotGrpDf = sqldf('SELECT *, group_concat(SYMBOL) AS SYMBOLS FROM annotDf GROUP BY Location')


## Expand gene sets and check associations
assColNames = c("distalDhsToPromoterDhs",
                "dhsToGeneExpression",
                "fantom5EnhancerTssAssociations",
                "allDbSuperEnhancerGeneAssociations",
                "fourDGenomeHomoSapiens", 
                "insitu_HiC_GM12878_100kb_MAPQGE30_100kb_SQRTVC",
                "insitu_HiC_HMEC_100kb_intra_MAPQGE30_100kb_SQRTVC",
                "insitu_HiC_HUVEC_100kb_intra_MAPQGE30_100kb_SQRTVC",
                "insitu_HiC_IMR90_100kb_intra_MAPQGE30_100kb_SQRTVC",
                "insitu_HiC_K562_100kb_intra_MAPQGE30_100kb_SQRTVC",
                "insitu_HiC_KBM7_100kb_intra_MAPQGE30_100kb_SQRTVC",
                "insitu_HiC_NHEK_100kb_intra_MAPQGE30_100kb_SQRTVC")

assGenesByLocLst = list()
for (i in 1:nrow(annotGrpDf)) {
    knownGenes = unique(strsplit(annotGrpDf[i,]$SYMBOLS, ",")[[1]], na.rm = T)
    knownGenes = knownGenes[knownGenes != ""]
    newGenes = c()
    for (assColName in assColNames) {
        if (!is.na(annotGrpDf[i, assColName])) {
            if (grepl("four", assColName)) {
                suppGenes = unique(unlist(sapply(sapply(strsplit(annotGrpDf[i, assColName], ",")[[1]],
                                                        function(x) strsplit(x, "_", fixed = T)[[1]][4]),
                                                 function(y) strsplit(y, "|", fixed = T)[[1]])))
                suppGenes = suppGenes[suppGenes != "NA"]
            } else if (grepl("insitu", assColName)) {
                suppGenes = unique(unlist(sapply(sapply(strsplit(annotGrpDf[i, assColName], ",")[[1]],
                                                        function(x) strsplit(x, "_", fixed = T)[[1]][7]),
                                                 function(y) strsplit(y, "|", fixed = T)[[1]])))
                suppGenes = suppGenes[suppGenes != "NA"]
            } else {
                suppGenes = unique(strsplit(annotGrpDf[i, assColName], ",")[[1]])
            }
            newGenes = c(newGenes, suppGenes)
        }
    }
    assGenesByLocLst[[annotGrpDf[i,]$Location]] = sort(unique(c(knownGenes, newGenes)))
}


## Test significance of expression changes of nearby genes and mutual exclusivity with other somatic mutations
assDf = data.frame()
for (i in 1:nrow(annotGrpDf)) {
    hotspotName = paste("chr", annotGrpDf[i, "Location"], sep = "")
    assGenes = assGenesByLocLst[[annotGrpDf[i, "Location"]]]
    if (length(assGenes) == 0) next
    cat(sprintf("Processing %d of %d (%d genes)...\n", i, nrow(annotGrpDf), length(assGenes)))
    assGrpDf <- foreach (assGene=assGenes, .combine = rbind) %dopar% {
        hotspotMtSids = as.vector(sapply(strsplit(annotGrpDf[i, "sids"], ",", fixed = T)[[1]],
                                         function(x) { 
                                             elems = strsplit(x, "-", fixed = T)[[1]][1:4]
                                             elems[1] = "TCGA"
                                             elems[4] = substring(elems[4], 1, 2)
                                             paste(elems, collapse = "_") }))
        hotspotWtSids = setdiff(totSids, hotspotMtSids)

        cnvByGeneGeneDf = cnvByGeneDf[assGene,]
        if (!all(is.na(cnvByGeneGeneDf)) & nrow(cnvByGeneGeneDf) > 0) {
            cnvMtSids = colnames(cnvByGeneGeneDf[, cnvByGeneGeneDf > 0.3 | cnvByGeneGeneDf < -0.3])
        } else {
            cnvMtSids = c()
        }

        ## Mutual exclusivity
        cnvMtSids = intersect(cnvMtSids, totSids)
        mafMtSids = intersect(unique(subset(mafPanFuncDf, Hugo_Symbol == assGene)$sid), totSids)
        somMtSids = union(mafMtSids, cnvMtSids)
        bothMutCnt = length(intersect(hotspotMtSids, somMtSids))
        hotspotOnlyMutCnt = length(setdiff(hotspotMtSids, somMtSids))
        somOnlyMutCnt = length(setdiff(somMtSids, hotspotMtSids))
        noneMutCnt = length(setdiff(hotspotWtSids, somMtSids))
        hotspotMat = matrix(c(bothMutCnt, somOnlyMutCnt, hotspotOnlyMutCnt, noneMutCnt),
                            nrow = 2, dimnames = list(Hotspot = c("MT", "WT"), Other = c("MT", "WT")))
        hotspotMutExcFisher = fisher.test(hotspotMat, alternative = "less")
        hotspotMutExcPvalue = hotspotMutExcFisher$p.value
        hotspotMutExcOr = as.numeric(hotspotMutExcFisher$estimate)

        if (hotspotMutExcPvalue < 0.05) {
            cat(sprintf("Found significant mutual exclusivity (odds ratio: %f, p-value: %f) between %s in %s mutations!\n",
                        hotspotMutExcOr, hotspotMutExcPvalue, hotspotName, assGene))
        }

        ## Expression changes
        geneExpLog2FoldGeneDt = geneExpLog2FoldPanDt[assGene,]
        if (nrow(geneExpLog2FoldGeneDt) > 1) {
            geneExpLog2FoldGeneDt = geneExpLog2FoldGeneDt[geneExpLog2FoldGeneDt$entrez == geneInfoExtDt[Symbol == assGene,]$GeneID[1],]
            if (nrow(geneExpLog2FoldGeneDt) > 1) {
                geneExpLog2FoldGeneDt = geneExpLog2FoldGeneDt[1, intersect(colnames(geneExpLog2FoldPanDt), totSids), with = F]
            }
        } else if (nrow(geneExpLog2FoldGeneDt) == 1) {
            geneExpLog2FoldGeneDt = geneExpLog2FoldGeneDt[, intersect(colnames(geneExpLog2FoldPanDt), totSids), with = F]
        }

        if (nrow(geneExpLog2FoldGeneDt) > 0) {
            expSids = grep("TCGA", names(geneExpLog2FoldGeneDt), value = T)
            expHotspotWtSids = setdiff(intersect(hotspotWtSids, expSids), cnvMtSids)
            expHotspotMtSids = intersect(hotspotMtSids, expSids)
            geneExpLog2FoldGeneMt = as.numeric(geneExpLog2FoldGeneDt[, expHotspotMtSids, with = F])
            geneExpLog2FoldGeneWt = as.numeric(geneExpLog2FoldGeneDt[, expHotspotWtSids, with = F])
            ttPvalueGeneLog2FoldExpMutVsWdt = my.t.test.p.value(geneExpLog2FoldGeneMt, geneExpLog2FoldGeneWt)
            wtPvalueGeneLog2FoldExpMutVsWdt = my.wilcox.test.p.value(geneExpLog2FoldGeneMt, geneExpLog2FoldGeneWt)
            if (!is.na(ttPvalueGeneLog2FoldExpMutVsWdt) & !is.na(wtPvalueGeneLog2FoldExpMutVsWdt)) {
                if (ttPvalueGeneLog2FoldExpMutVsWdt < 0.05 | wtPvalueGeneLog2FoldExpMutVsWdt < 0.05) {
                    cat(sprintf("Found significant expression change of %s in %s mutant group!\n", assGene, hotspotName))
                    geneLog2FoldExpMutVsWdtDf = data.frame()
                    geneLog2FoldExpMutVsWdtDf = rbind(geneLog2FoldExpMutVsWdtDf, data.frame(type=rep("MT", length(expHotspotMtSids)), log2FoldRsem = geneExpLog2FoldGeneMt))
                    geneLog2FoldExpMutVsWdtDf = rbind(geneLog2FoldExpMutVsWdtDf, data.frame(type=rep("WT", length(expHotspotWtSids)), log2FoldRsem = geneExpLog2FoldGeneWt))
                    p = ggplot(geneLog2FoldExpMutVsWdtDf, aes(factor(type), log2FoldRsem)) +
                        geom_boxplot() + 
                        theme(plot.title = element_text(size = baseFontSize + 2, face="bold"),
                            axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
                            axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
                            axis.text.x = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 0, hjust = NULL),
                            axis.text.y = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_rect(color = "black", fill="white"),
                            legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
                            legend.text = element_text(size = baseFontSize, family="sans"),
                            legend.direction = "horizontal",
                            legend.position = "bottom") +
                        scale_x_discrete(name=sprintf("\nWilcoxon signed-rank test: %f\nStudent's t-test: %f",
                                                      wtPvalueGeneLog2FoldExpMutVsWdt, ttPvalueGeneLog2FoldExpMutVsWdt)) +
                        scale_y_continuous(name = sprintf("Log2 fold change of %s expression\n", assGene))
                    mutVsWdtExpPlotDir = file.path(figDir, "expression", annotGrpDf[i, "space"])
                    dir.create(mutVsWdtExpPlotDir, recursive = TRUE, showWarnings = FALSE)
                    mutVsWdtExpPlotFile = file.path(mutVsWdtExpPlotDir, sprintf("%s-%s-MutVsWdt-Log2FoldRsem.pdf", hotspotName, assGene))
                    ggsave(filename = mutVsWdtExpPlotFile, plot = p, width = 5, height = 6)
                }
            }
        } else {
            ttPvalueGeneLog2FoldExpMutVsWdt = NA
            wtPvalueGeneLog2FoldExpMutVsWdt = NA
        }
        data.frame(Location = annotGrpDf[i,]$Location,
                   SYMBOL = assGene,
                   mutualExclusivityPvalue = hotspotMutExcPvalue,
                   mutualExclusivityLog2OddsRatio = log2(hotspotMutExcOr),
                   log2FoldExpressionTtestPvalue = ttPvalueGeneLog2FoldExpMutVsWdt,
                   log2FoldExpressionWilcoxPvalue = wtPvalueGeneLog2FoldExpMutVsWdt)
    }
    assDf = rbind(assDf, assGrpDf)
}


## Merge expanded gene associations with original annotation results
annotGrpDf$SYMBOL = NULL
annotGrpDf$SYMBOLS = NULL
assAnnotDf = merge(assDf, annotGrpDf, by = c("Location"), all.x = T)

## Check which method supports the associations
for (i in 1:nrow(assAnnotDf)) {
    gene = assAnnotDf[i,]$SYMBOL
    if (gene != "") {
        suppColNames = c()
        for (assColName in assColNames) {
            if (grepl(gene, assAnnotDf[i, assColName], fixed = T)) {
                cat(sprintf("Association between hotspot (%s) and %s is supported by %s\n",
                            assAnnotDf[i, "Location"], gene, assColName))
                suppColNames = c(suppColNames, assColName)
            }
        }
        assAnnotDf[i, "supported_by"] = paste(suppColNames, collapse = ";")
    }
}
assAnnotDf[is.na(assAnnotDf$supported_by), "supported_by"] = ""
for (assColName in assColNames) { assAnnotDf[, assColName] = NULL }


## Save the results!
assFile = gsub("txt", "asstested.txt", annotFile, fixed = T)
write.table(assAnnotDf, assFile, row.names = F, col.names = T, sep = "\t", quote = F)

