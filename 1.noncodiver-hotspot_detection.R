args = commandArgs(TRUE)
scntFile = args[1]
scntSelGrpFile = args[2]
scntSelGrpVepFile = args[3]
distCut = as.integer(args[4])
hotspotMargin = as.integer(args[5])
poisMargin = as.integer(args[6])
numCores = as.integer(args[7])

#scntFile = "/n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/hotspot/TCGA_16_Cancer_Types.wgs.somatic.chr5.sanitized.scnt.txt"
#scntSelGrpFile = "/n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/hotspot/TCGA_16_Cancer_Types.wgs.somatic.chr5.sanitized.scnt.hotspot100.txt"
#scntSelGrpVepFile = "/n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/hotspot/TCGA_16_Cancer_Types.wgs.somatic.chr5.sanitized.scnt.hotspot100.vep_in.txt"
#distCut = 100
#poisMargin = 10000
#numCores = 10

require(doMC)
registerDoMC(numCores)

require(gdata)
require(gtools)
require(Biostrings)
require(GenomicRanges)
require(BSgenome.Hsapiens.UCSC.hg19)


##
## Initialize global variables
##
rootDir = "/home/sl279"
baseDir = file.path(rootDir, "BiO/Research/NoncoDiver")
vcfDir = file.path(baseDir, "vcf")
panVcfDir = file.path(vcfDir, "pancan")
figDir = file.path(baseDir, "figure")

chrs = c(1:22, "X", "Y")
chrs = factor(chrs, level=chrs)
cchrs = sapply(chrs, function(x) paste("chr", x, sep=""))
cchrs = factor(cchrs, level=cchrs)

# Load chrom size information
hg19File = file.path(baseDir, "ucsc/database/hg19.genome")
hg19Df = read.delim(hg19File, header=T, as.is=T)
hg19Df = hg19Df[hg19Df$chrom %in% chrs,]
hg19Df$chrom = factor(hg19Df$chrom, levels=chrs)
hg19Df = hg19Df[order(hg19Df$chrom),]

# Cluster mutations
cat(sprintf("Processing %s ...\n", scntFile))
scntDf = read.delim(scntFile, header = F, as.is = T)
colnames(scntDf) = c("chrom", "pos", "ref", "alt", "scnt", "sids", "cadd_raw", "cadd_phred")

scntGrp = list()
grpIdx = 1
for (i in 1:nrow(scntDf)) {
    if (i == 1) {
        scntGrp[[grpIdx]] = i
    } else if (scntDf[i, "pos"] - scntDf[i-1, "pos"] <= distCut) {
        scntGrp[[grpIdx]] = c(scntGrp[[grpIdx]], i)
    } else {
        grpIdx = grpIdx + 1
        scntGrp[[grpIdx]] = i
    }
}

scntGrpSids = mclapply(scntGrp, function(x) { unique(unlist(sapply(scntDf[x, "sids"], function(x) { unlist(strsplit(strsplit(x, ",")[[1]], ";")) }))) }, mc.cores = numCores)
scntGrpSids = mclapply(scntGrpSids, function(x) { paste(x, collapse = ",") }, mc.cores = numCores)
scntGrpSids = unlist(scntGrpSids)
scntGrpSidCnt = as.vector(sapply(scntGrpSids, function(x) length(strsplit(x, ",", fixed = T)[[1]])))
scntSelGrp = scntGrp[scntGrpSidCnt > 1]
scntSelGrpSids = scntGrpSids[scntGrpSidCnt > 1]
scntSelGrpSidCnt = scntGrpSidCnt[scntGrpSidCnt > 1]
scntSelGrpStartPos = sapply(scntSelGrp, function(x) { min(scntDf[x, "pos"]) - hotspotMargin })
scntSelGrpEndPos = sapply(scntSelGrp, function(x) { max(scntDf[x, "pos"]) + hotspotMargin })
scntSelGrpWidth = scntSelGrpEndPos - scntSelGrpStartPos + 1
scntSelGrpNumSnvLoci = sapply(scntSelGrp, function(x) length(x))
scntSelGrpSnvCnt = sapply(scntSelGrp, function(x) { sum( as.numeric(unlist(sapply(scntDf[x, "scnt"], function(x) strsplit(x, ",")[[1]])))) } )
chrL = hg19Df[hg19Df$chrom == scntDf$chrom[1],]$size
scntSelGrpPvalues = mclapply(1:length(scntSelGrp),
                            function(i) {
                                grpB = scntSelGrpStartPos[i]
                                grpE = scntSelGrpEndPos[i]
                                grpL = grpE - grpB + 1
                                grpS = scntSelGrpSidCnt[i]
                                smpB = grpB - poisMargin
                                smpE = grpE + poisMargin
                                smpB = ifelse(smpB < 0, 0, smpB)
                                smpE = ifelse(smpE > chrL, chrL, smpE)
                                smpL = smpE - smpB + 1
                                smpS = length(unique(unlist(sapply(unlist(sapply(subset(scntDf, pos >= smpB & pos <= smpE)$sids,
                                                                                 function(x) { strsplit(x, ",")[[1]] })),
                                                                   function(y) { strsplit(y, ";")[[1]] }))))
                                lambda = smpS / smpL
                                grpPoiP = poisson.test(grpS, lambda * grpL, alternative = "greater")$p.value
                                return(grpPoiP)
                            }, mc.cores = numCores)
scntSelGrpPvalues = unlist(scntSelGrpPvalues)
scntSelGrpChrFdrs = p.adjust(scntSelGrpPvalues, method = "fdr")
scntSelGrpNrSidRef = getSeq(BSgenome.Hsapiens.UCSC.hg19,
                        GRanges(seqnames = Rle(gsub("(.*)", "chr\\1", rep(scntDf$chrom[1], length(scntSelGrp)), perl = T)),
                                ranges = IRanges(start = scntSelGrpStartPos, end = scntSelGrpEndPos),
                                strand = "+"),
                        as.character = T)
scntSelGrpNrSidAlt = rep("-", length(scntSelGrp))
scntSelGrpCaddRaw = sapply(scntSelGrp, function(x) { mean(as.numeric(unlist(sapply(scntDf[x, "cadd_raw"], function(y) strsplit(y, ",", fixed = T)[[1]])))) } )
scntSelGrpCaddPhred = sapply(scntSelGrp, function(x) { mean(as.numeric(unlist(sapply(scntDf[x, "cadd_phred"], function(y) strsplit(y, ",", fixed = T)[[1]])))) } )
scntSelGrpDf = data.frame(space = rep(scntDf[1, "chrom"], length(scntSelGrp)),
                          start = scntSelGrpStartPos,
                          end = scntSelGrpEndPos,
                          width = scntSelGrpWidth,
                          mloci = scntSelGrpNumSnvLoci,
                          ref = scntSelGrpNrSidRef,
                          alt = scntSelGrpNrSidAlt,
                          allele = paste(scntSelGrpNrSidRef, scntSelGrpNrSidAlt, sep = "/"),
                          strand = "+",
                          mcnt = scntSelGrpSnvCnt,
                          mcnt_per_hsize = scntSelGrpSnvCnt / scntSelGrpWidth,
                          mcnt_per_mloci = scntSelGrpSnvCnt / scntSelGrpNumSnvLoci,
                          scnt = scntSelGrpSidCnt,
                          scnt_p_value = scntSelGrpPvalues,
                          chr_adj_scnt_p_value = scntSelGrpChrFdrs,
                          avg_cadd_raw = scntSelGrpCaddRaw,
                          avg_cadd_phred = scntSelGrpCaddPhred,
                          sids = scntSelGrpSids)
scntSelGrpDf = scntSelGrpDf[order(scntSelGrpDf$chr_adj_scnt_p_value),]

write.table(scntSelGrpDf, scntSelGrpFile, row.names = F, col.names = T, sep = "\t", quote = F)
write.table(scntSelGrpDf[, c("space", "start", "end", "allele", "strand")],
            scntSelGrpVepFile, row.names = F, col.names = F, sep = "\t", quote = F)

