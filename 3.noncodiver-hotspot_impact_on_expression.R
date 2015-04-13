args = commandArgs(TRUE)
annotFile = args[1]
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
require(data.table)
require(Biostrings)
require(GenomicRanges)
require(org.Hs.eg.db)
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

## Load mRNA expression data from GDAC
rsemLog2FoldPanRdata = file.path(gdacStdDataDir, "rsemLog2FoldPan.RData")
load(rsemLog2FoldPanRdata)

## Load SCNA data for filter
cnvByGeneRdata = file.path(gdacStdDataDir, "cnvByGene.RData")
load(cnvByGeneRdata)

## Test significance of expression changes of nearby genes
annotDf = read.delim(annotFile, header = T, as.is = T)
expTestedDf <- foreach (i=1:nrow(annotDf), .combine=rbind) %dopar% {
    hotspotName = sprintf("%s_%s_%s",
                          annotDf[i, "space"],
                          annotDf[i, "start"], 
                          annotDf[i, "end"])
    cat(sprintf("Processing %d of %d ...\n", i, nrow(annotDf)))
    geneName = annotDf[i, "SYMBOL"]
    cnvByGeneGeneDf = cnvByGeneDf[geneName,]
    if (!all(is.na(cnvByGeneGeneDf)) & nrow(cnvByGeneGeneDf) > 0) {
        cpids = colnames(cnvByGeneGeneDf[, cnvByGeneGeneDf > 0.3 | cnvByGeneGeneDf < -0.3])
    } else {
        cpids = c()
    }
    rsemLog2FoldGeneDf = subset(rsemLog2FoldPanDf, symbol == geneName)
    if (nrow(rsemLog2FoldGeneDf) > 0) {
        mpids = as.vector(sapply(strsplit(annotDf[i, "sids"], ",", fixed = T)[[1]],
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
        annotDf[i, "ttPvalueGeneLog2FoldExpMutVsWdt"] = ttPvalueGeneLog2FoldExpMutVsWdt
        annotDf[i, "wtPvalueGeneLog2FoldExpMutVsWdt"] = wtPvalueGeneLog2FoldExpMutVsWdt
        if (!is.na(ttPvalueGeneLog2FoldExpMutVsWdt) & !is.na(wtPvalueGeneLog2FoldExpMutVsWdt)) {
            if (ttPvalueGeneLog2FoldExpMutVsWdt <= 0.05 | wtPvalueGeneLog2FoldExpMutVsWdt <= 0.05) {
                cat(sprintf("Found a significant expression change of %s in %s mutant group!\n", geneName, hotspotName))
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
                        legend.position = "bottom"
                        ) +
                    scale_x_discrete(name=sprintf("\nWilcoxon signed-rank test: %f\nStudent's t-test: %f", wtPvalueGeneLog2FoldExpMutVsWdt, ttPvalueGeneLog2FoldExpMutVsWdt)) +
                    scale_y_continuous(name = sprintf("Log2 fold change of %s expression\n", geneName))
                mutVsWdtExpPlotFile = file.path(figDir, sprintf("%s-%s-MutVsWdt-Log2FoldRsem.pdf", hotspotName, geneName))
                ggsave(filename = mutVsWdtExpPlotFile, plot = p, width = 5, height = 6)
            }
        }
    } else {
        annotDf[i, "ttPvalueGeneLog2FoldExpMutVsWdt"] = NA
        annotDf[i, "wtPvalueGeneLog2FoldExpMutVsWdt"] = NA
    }
    annotDf[i,]
}

expTestedFile = gsub("annotated.txt", "annotated.exptested.txt", annotFile, fixed = T)
write.table(expTestedDf, expTestedFile, row.names = F, col.names = T, sep = "\t", quote = F)

