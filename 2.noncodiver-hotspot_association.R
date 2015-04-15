args = commandArgs(TRUE)
annotFile = args[1]
annotFile = "/n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/hotspot/TCGA_16_Cancer_Types.wgs.somatic.chr10.sanitized.scnt.hotspot100.fdr0.05.annotated.txt"
numCores = 10

##
## Load libraries
##
require(doMC)
registerDoMC(numCores)

require(proxy)
require(gdata)
require(gtools)
require(ggplot2)
require(annotate)
require(snpStats)
require(Biostrings)
require(org.Hs.eg.db)
require(GenomicRanges)
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


# Load exome SNV/indel data from GDAC
mafRdata = file.path(gdacAnlDataDir, "mafPanFunc.RData")
load(mafRdata)

## Load mRNA expression data from GDAC
rsemLog2FoldPanRdata = file.path(gdacStdDataDir, "rsemLog2FoldPan.RData")
load(rsemLog2FoldPanRdata)

## Load SCNA data for filter
cnvByGeneRdata = file.path(gdacStdDataDir, "cnvByGene.RData")
load(cnvByGeneRdata)

## Load hotspot data
annotDf = read.delim(annotFile, header = T, as.is = T)
#nrow(annotDf)
#names(annotDf)

annotDf$allele = NULL
annotDf$strand = NULL
annotGrpDf = sqldf('SELECT *, group_concat(SYMBOL) AS SYMBOLS FROM annotDf GROUP BY Location')
#nrow(annotGrpDf)
#head(annotGrpDf)

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
    newGenes = unique(newGenes)
    assGenes = unique(c(knownGenes, newGenes))
    assGenesByLocLst[[annotGrpDf[i,]$Location]] = assGenes
}

for (i in 1:nrow(annotDf)) {
    gene = annotDf[i,]$SYMBOL
    if (gene == "") {
        next
    } else {
        suppColNames = c()
        for (assColName in assColNames) {
            if (grepl(gene, annotDf[i, assColName], fixed = T)) {
                cat(sprintf("Association between hotspot (%s) and %s is supported by %s\n",
                            annotDf[i, "Location"], gene, assColName))
                suppColNames = c(suppColNames, assColName)
            }
        }
        annotDf[i, "supported_by"] = paste(suppColNames, collapse = ";")
    }
}
annotDf[is.na(annotDf$supported_by), "supported_by"] = ""


## Test significance of expression changes of nearby genes and mutual exclusivity with other somatic mutations
assDf = data.frame()
for (i = 1:nrow(annotGrpDf)) {
    hotspotName = paste("chr", annotGrpDf[i, "Location"], sep = "")
    cat(sprintf("Processing %d of %d ...\n", i, nrow(annotGrpDf)))
    assGenes = assGenesByLocLst[[annotGrpDf[i, "Location"]]]

    assGrpDf <- foreach (assGene=assGenes) %dopar% {
        cnvByGeneGeneDf = cnvByGeneDf[assGene,]
        if (!all(is.na(cnvByGeneGeneDf)) & nrow(cnvByGeneGeneDf) > 0) {
            cpids = colnames(cnvByGeneGeneDf[, cnvByGeneGeneDf > 0.3 | cnvByGeneGeneDf < -0.3])
        } else {
            cpids = c()
        }
        rsemLog2FoldGeneDf = subset(rsemLog2FoldPanDf, symbol == assGene)
        if (nrow(rsemLog2FoldGeneDf) > 0) {
            mpids = as.vector(sapply(strsplit(annotGrpDf[i, "sids"], ",", fixed = T)[[1]],
                                    function(x) { 
                                        elems = strsplit(x, "-", fixed = T)[[1]][2:4]
                                        elems[3] = substring(elems[3], 1, 2)
                                        paste(elems, collapse = "__")
                                    }))
            epids = grep("__", names(rsemLog2FoldPanDf), value = T)
            wepids = setdiff(setdiff(epids, mpids), cpids)
            mepids = intersect(mpids, epids)
            rsemLog2FoldGeneMt = as.numeric(as.vector(rsemLog2FoldGeneDf[, mepids]))
            rsemLog2FoldGeneWt = as.numeric(as.vector(rsemLog2FoldGeneDf[, wepids]))
            ttPvalueGeneLog2FoldExpMutVsWdt = my.t.test.p.value(rsemLog2FoldGeneMt, rsemLog2FoldGeneWt)
            wtPvalueGeneLog2FoldExpMutVsWdt = my.wilcox.test.p.value(rsemLog2FoldGeneMt, rsemLog2FoldGeneWt)
            if (!is.na(ttPvalueGeneLog2FoldExpMutVsWdt) & !is.na(wtPvalueGeneLog2FoldExpMutVsWdt)) {
                if (ttPvalueGeneLog2FoldExpMutVsWdt <= 0.05 | wtPvalueGeneLog2FoldExpMutVsWdt <= 0.05) {
                    cat(sprintf("Found a significant expression change of %s in %s mutant group!\n", assGene, hotspotName))
                    geneLog2FoldExpMutVsWdtDf = data.frame()
                    geneLog2FoldExpMutVsWdtDf = rbind(geneLog2FoldExpMutVsWdtDf, data.frame(type=rep("MT", length(mepids)), log2FoldRsem = rsemLog2FoldGeneMt))
                    geneLog2FoldExpMutVsWdtDf = rbind(geneLog2FoldExpMutVsWdtDf, data.frame(type=rep("WT", length(wepids)), log2FoldRsem = rsemLog2FoldGeneWt))
                    p = ggplot(geneLog2FoldExpMutVsWdtDf, aes(factor(type), log2FoldRsem)) +
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
                        scale_x_discrete(name=sprintf("\nWilcoxon signed-rank test: %f\nStudent's t-test: %f", wtPvalueGeneLog2FoldExpMutVsWdt, ttPvalueGeneLog2FoldExpMutVsWdt)) +
                        scale_y_continuous(name = sprintf("Log2 fold change of %s expression\n", assGene))
                    mutVsWdtExpPlotFile = file.path(figDir, sprintf("%s-%s-MutVsWdt-Log2FoldRsem.pdf", hotspotName, assGene))
                    ggsave(filename = mutVsWdtExpPlotFile, plot = p, width = 5, height = 6)
                }
            }
        } else {
            ttPvalueGeneLog2FoldExpMutVsWdt = NA
            wtPvalueGeneLog2FoldExpMutVsWdt = NA
        }

        return(data.frame(Location = annotGrpDf[i,]$Location,
                          SYMBOL = assGene,
                          ttPvalueGeneLog2FoldExpMutVsWdt = ttPvalueGeneLog2FoldExpMutVsWdt,
                          wtPvalueGeneLog2FoldExpMutVsWdt = wtPvalueGeneLog2FoldExpMutVsWdt))
    }
    assDf = rbind(assDf, assGrpDf)
}

assFile = gsub("annotated.txt", "annotated.asstested.txt", annotFile, fixed = T)
write.table(expTestedDf, expTestedFile, row.names = F, col.names = T, sep = "\t", quote = F)

