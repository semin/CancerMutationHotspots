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

refGeneFile = file.path("/home/sl279/BiO/Research/NoncoDiver/ucsc/database/refGene.txt.gz")
refGeneDf = read.delim(gzfile(refGeneFile), header = F, as.is = T)[, c(3:6,13)]
colnames(refGeneDf) = c("space", "strand", "start", "end", "symbol")
refGeneDf = refGeneDf[refGeneDf$space %in% cchrs,]
refGeneDf = with(refGeneDf, refGeneDf[order(space, start, end),])
refGeneGr = with(refGeneDf, GRanges(seqnames = Rle(space), IRanges(start = start, end = end), strand = strand, symbol = symbol))

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

targetGrA = sort(GRanges(seqnames = Rle(ichrs[1]),
                   IRanges(start = unique(obsDf$i + 1), end = unique(obsDf$i + binSize)),
                   strand = "*"))
geneHitsA = findOverlaps(targetGrA, refGeneGr)
for(i in unique(queryHits(geneHitsA))) {
    values(targetGrA)[i, "symbols"] = paste(unique(refGeneGr[subjectHits(geneHitsA[queryHits(geneHitsA) == i])]$symbol), collapse = "|")
}
targetDfA = as.data.frame(targetGrA)

targetGrB = sort(GRanges(seqnames = Rle(ichrs[length(ichrs)]),
                   IRanges(start = unique(obsDf$j + 1), end = unique(obsDf$j + binSize)),
                   strand = "*"))
geneHitsB = findOverlaps(targetGrB, refGeneGr)
for(i in unique(queryHits(geneHitsB))) {
    values(targetGrB)[i, "symbols"] = paste(unique(refGeneGr[subjectHits(geneHitsB[queryHits(geneHitsB) == i])]$symbol), collapse = "|")
}
targetDfB = as.data.frame(targetGrB)

obsDf$chrA = ichrs[1]
obsDf$startA = as.integer(obsDf$i + 1)
obsDf$endA = as.integer(obsDf$i + binSize)
obsDf$chrB = ichrs[length(ichrs)]
obsDf$startB = as.integer(obsDf$j + 1)
obsDf$endB = as.integer(obsDf$j + binSize)
obsDf = sqldf('SELECT o.chrA, o.startA, o.endA, o.chrB, o.startB, o.endB, o.c, o.nc, o.loe, t.symbols AS symbolsA
              FROM obsDf AS o
              LEFT JOIN targetDfA AS t
              ON o.chrA = t.seqnames AND o.startA = t.start AND o.endA = t.end')
obsDf = sqldf('SELECT o.chrA, o.startA, o.endA, symbolsA, o.chrB, o.startB, o.endB, o.c, o.nc, o.loe, t.symbols AS symbolsB
              FROM obsDf AS o
              LEFT JOIN targetDfB AS t
              ON o.chrB = t.seqnames AND o.startB = t.start AND o.endB = t.end')
#obsSigDf = obsDf[obsDf$nc > as.numeric(quantile(obsDf$nc, 0.95, na.rm = T)) &
                 #obsDf$c > as.numeric(quantile(obsDf$c, 0.95, na.rm = T)),]
obsSigDf = obsDf[obsDf$c > as.numeric(quantile(obsDf$c, 0.95)),]
bedDf = with(obsSigDf, rbind(data.frame(space = chrA, start = startA, end = endA, partner = paste(chrB, startB, endB, c, nc, loe, symbolsB, sep = "_")),
                             data.frame(space = chrB, start = startB, end = endB, partner = paste(chrA, startA, endA, c, nc, loe, symbolsA, sep = "_"))))
bedGrpDf = sqldf('SELECT * GROUP BY space, start, end ORDER BY space, start, end ASC')
write.table(bedGrpDf, outFile, sep = "\t", row.names=F, col.names=F, quote=FALSE)

