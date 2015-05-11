#binSize = 100000
#obsFile = "/home/sl279/BiO/Research/Hotspot/insitu-hic/GM12878_combined/100kb_resolution_interchromosomal/chr7_chr17/MAPQGE30/chr7_17_100kb.RAWobserved"
#nrm1File = "/home/sl279/BiO/Research/Hotspot/insitu-hic/GM12878_combined/100kb_resolution_interchromosomal/chr7_chr17/MAPQGE30/chr7_100kb.SQRTVCnorm"
#nrm2File = "/home/sl279/BiO/Research/Hotspot/insitu-hic/GM12878_combined/100kb_resolution_interchromosomal/chr7_chr17/MAPQGE30/chr17_100kb.SQRTVCnorm"
#expFile = "NA"
#outFile = "/home/sl279/BiO/Research/Hotspot/insitu-hic/GM12878_combined/100kb_resolution_interchromosomal/chr7_chr17/MAPQGE30/chr7_17_100kb.SQRTVC.bed"

args = commandArgs(TRUE)
binSize = as.numeric(args[1])
obsFile = args[2]
nrm1File = args[3]
nrm2File = args[4]
expFile = args[5]
outFile = args[6]

require(sqldf)
require(GenomicRanges)

chrs = c(1:22, "X", "Y")
chrs = factor(chrs, level=chrs)
cchrs = sapply(chrs, function(x) paste("chr", x, sep=""))
cchrs = factor(cchrs, level=cchrs)

refGeneFile = file.path("/home/sl279/BiO/Research/Hotspot/ucsc/database/refGene.txt.gz")
refGeneDf = read.delim(gzfile(refGeneFile), header = F, as.is = T)[, c(2:6,13)]
colnames(refGeneDf) = c("transcript", "space", "strand", "start", "end", "gene")
refGeneDf$tss = with(refGeneDf, ifelse(strand == "+", start, end))
refGenePromDf = sqldf('SELECT gene, strand, space, tss FROM refGeneDf GROUP BY gene, tss')
prom_margin = 2000
refGenePromDf$prom_start = refGenePromDf$tss - prom_margin
refGenePromDf$prom_end = refGenePromDf$tss + prom_margin
refGenePromDf$space = factor(refGenePromDf$space, levels = cchrs)
refGenePromDf = refGenePromDf[refGenePromDf$space %in% cchrs,]
refGenePromGr = with(refGenePromDf, sort(GRanges(seqnames = Rle(space),
                                            IRanges(start = prom_start, end = prom_end),
                                            strand = strand,
                                            gene = gene,
                                            tss = tss)))

obsDf = read.delim(obsFile, header = F, col.names = c("i", "j", "c"))
nrm1Df = read.delim(nrm1File, header = F, col.names = c("n"))

if (basename(nrm2File) != "NA") {
    nrm2Df = read.delim(nrm2File, header = F, col.names = c("n"))
    obsDf$nc = obsDf$c / (nrm1Df[obsDf$i / binSize + 1,] * nrm2Df[obsDf$j / binSize + 1,])
} else {
    obsDf$nc = obsDf$c / (nrm1Df[obsDf$i / binSize + 1,] * nrm1Df[obsDf$j / binSize + 1,])
}

if (basename(expFile) != "NA") {
    expDf = read.delim(expFile, header = F, col.names = c("e"))
    obsDf$oe = obsDf$nc / expDf[(abs(obsDf$j - obsDf$i) / binSize + 1),]
    obsDf$loe = log2(obsDf$oe)
} else {
    obsDf$oe = NA
    obsDf$loe = NA
}

ichrs = strsplit(basename(dirname(dirname(obsFile))), "_", fixed = T)[[1]]

targetGrA = unique(sort(GRanges(seqnames = Rle(ichrs[1]),
                                IRanges(start = unique(obsDf$i + 1), end = unique(obsDf$i + binSize)),
                                strand = "*")))
geneHitsA = findOverlaps(targetGrA, refGenePromGr)
for(i in unique(queryHits(geneHitsA))) {
    values(targetGrA)[i, "genes"] = paste(unique(refGenePromGr[subjectHits(geneHitsA[queryHits(geneHitsA) == i])]$gene), collapse = "|")
}
targetDfA = as.data.frame(targetGrA)

targetGrB = unique(sort(GRanges(seqnames = Rle(ichrs[length(ichrs)]),
                                IRanges(start = unique(obsDf$j + 1), end = unique(obsDf$j + binSize)),
                                strand = "*")))
geneHitsB = findOverlaps(targetGrB, refGenePromGr)
for(i in unique(queryHits(geneHitsB))) {
    values(targetGrB)[i, "genes"] = paste(unique(refGenePromGr[subjectHits(geneHitsB[queryHits(geneHitsB) == i])]$gene), collapse = "|")
}
targetDfB = as.data.frame(targetGrB)

obsDf$chrA = ichrs[1]
obsDf$startA = as.integer(obsDf$i + 1)
obsDf$endA = as.integer(obsDf$i + binSize)
obsDf$chrB = ichrs[length(ichrs)]
obsDf$startB = as.integer(obsDf$j + 1)
obsDf$endB = as.integer(obsDf$j + binSize)
obsDf = sqldf('SELECT o.chrA, o.startA, o.endA, o.chrB, o.startB, o.endB, o.c, o.nc, o.loe, t.genes AS genesA
              FROM obsDf AS o
              LEFT JOIN targetDfA AS t
              ON o.chrA = t.seqnames AND o.startA = t.start AND o.endA = t.end')
obsDf = sqldf('SELECT o.chrA, o.startA, o.endA, genesA, o.chrB, o.startB, o.endB, o.c, o.nc, o.loe, t.genes AS genesB
              FROM obsDf AS o
              LEFT JOIN targetDfB AS t
              ON o.chrB = t.seqnames AND o.startB = t.start AND o.endB = t.end')
obsSigDf = obsDf[obsDf$nc > as.numeric(quantile(obsDf$nc, 0.999)) & obsDf$c > as.numeric(quantile(obsDf$c, 0.999)),]
bedDf = with(obsSigDf, rbind(data.frame(space = chrA, start = startA, end = endA, partner = paste(chrB, startB, endB, c, nc, loe, genesB, sep = "_")),
                             data.frame(space = chrB, start = startB, end = endB, partner = paste(chrA, startA, endA, c, nc, loe, genesA, sep = "_"))))
bedDf = with(bedDf, bedDf[order(space, start, end),])
write.table(bedDf, outFile, sep = "\t", row.names=F, col.names=F, quote=FALSE)

