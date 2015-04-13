rm(list=ls())

require(Hmisc)
require(GenomicRanges)

baseDir = file.path("/home/sl279/BiO/Research/NoncoDiver")

chrs = c(1:22, "X", "Y")
chrs = factor(chrs, level=chrs)
cchrs = sapply(chrs, function(x) paste("chr", x, sep=""))
cchrs = factor(cchrs, level=cchrs)

##
## Merge regulatory dnase bed files
##

for (type in c("dyadic", "enh", "prom")) {
    print(type)
    totBedDf = data.frame()
    bedFiles = Sys.glob(file.path(baseDir, sprintf("roadmap/dnase/BED_files_%s/*.bed.gz", type)))
    for (bedFile in bedFiles) {
        print(bedFile)
        eid = gsub(".*(E\\d+)\\.bed\\.gz", "\\1", basename(bedFile), perl = T)
        bedDf = read.delim(gzfile(bedFile), header = F, as.is = T)[, c(1:3)]
        colnames(bedDf) = c("space", "start", "end")
        bedDf$eid = eid
        totBedDf = rbind(totBedDf, bedDf)
    }

    totBedDf$space = factor(totBedDf$space, levels = cchrs)
    totBedDf = totBedDf[order(totBedDf$space, totBedDf$start, totBedDf$end),]
    totBedGr = GRanges(seqnames = Rle(totBedDf$space),
                    IRanges(start = totBedDf$start,
                            end = totBedDf$end),
                    strand = "*",
                    eid = totBedDf$eid)

    totBedReducedGr = reduce(totBedGr)
    totBedReducedGr$ecnt = countOverlaps(totBedReducedGr, totBedGr)
    #ol = findOverlaps(totBedReducedGr, totBedGr)
    #totBedReducedGr$eids = unlist(mclapply(unique(queryHits(ol)),
                                        #function (i) {
                                            #paste(totBedGr[subjectHits(ol[queryHits(ol) == i])]$eid, collapse = ",")
                                        #}, mc.cores = 10))
    ctype = capitalize(type)
    #totBedReducedDf = as.data.frame(totBedReducedGr)[, c(1,2,3,6,7)]
    totBedReducedDf = as.data.frame(totBedReducedGr)[, c(1,2,3,6)]
    totBedReducedFile = file.path(baseDir, sprintf("roadmap/dnase/roadmapDnase%s.bed", ctype))
    write.table(totBedReducedDf, totBedReducedFile, sep = "\t", row.names=F, col.names=F, quote=FALSE)
}


##
## Merge peaks
##

peakFiles = Sys.glob(file.path(baseDir, "roadmap/peaks/gappedPeak/*.gappedPeak.gz"))
peakTypes = unique(sapply(peakFiles, function(x) { gsub("\\S+-(\\S+?)\\.*gappedPeak\\.gz", "\\1", basename(x)) }))

for (peakType in peakTypes) {
    print(peakType)
    markPeakTotDf = data.frame()
    markPeakFiles = Sys.glob(file.path(baseDir, sprintf("roadmap/peaks/gappedPeak/*%s.gappedPeak.gz", peakType)))
    for (markPeakFile in markPeakFiles) {
        print(markPeakFile)
        eid = gsub("(E\\d+)-.*", "\\1", basename(markPeakFile), perl = T)
        markPeakDf = read.delim(gzfile(markPeakFile), header = F, as.is = T)[, c(1:3)]
        colnames(markPeakDf) = c("space", "start", "end")
        markPeakDf = markPeakDf[markPeakDf$space %in% cchrs,]
        markPeakDf$eid = eid
        markPeakTotDf = rbind(markPeakTotDf, markPeakDf)
    }
    markPeakTotDf$space = factor(markPeakTotDf$space, levels = cchrs)
    markPeakTotDf = markPeakTotDf[order(markPeakTotDf$space, markPeakTotDf$start, markPeakTotDf$end),]
    markPeakTotGr = GRanges(seqnames = Rle(markPeakTotDf$space),
                            IRanges(start = markPeakTotDf$start,
                                    end = markPeakTotDf$end),
                            strand = "*",
                            eid = markPeakTotDf$eid)

    markPeakTotReducedGr = reduce(markPeakTotGr)
    markPeakTotReducedGr$ecnt = countOverlaps(markPeakTotReducedGr, markPeakTotGr)
    #ol = findOverlaps(markPeakTotReducedGr, markPeakTotGr)
    #markPeakTotReducedGr$eids = unlist(mclapply(unique(queryHits(ol)),
                                        #function (i) {
                                            #paste(markPeakTotGr[subjectHits(ol[queryHits(ol) == i])]$eid, collapse = ",")
                                        #}, mc.cores = 10))

    cPeakType = capitalize(peakType)
    #markPeakTotReducedDf = as.data.frame(markPeakTotReducedGr)[, c(1,2,3,6,7)]
    markPeakTotReducedDf = as.data.frame(markPeakTotReducedGr)[, c(1,2,3,6)]
    markPeakTotReducedFile = file.path(baseDir, sprintf("roadmap/peaks/roadmap%sGappedPeaks.bed", cPeakType))
    write.table(markPeakTotReducedDf, markPeakTotReducedFile, sep = "\t", row.names=F, col.names=F, quote=FALSE)
}


peakFiles = Sys.glob(file.path(baseDir, "roadmap/peaks/narrowPeak/*.narrowPeak.gz"))
peakTypes = unique(sapply(peakFiles, function(x) { gsub("\\S+-(\\S+?)\\.*narrowPeak\\.gz", "\\1", basename(x)) }))

for (peakType in peakTypes) {
    print(peakType)
    markPeakTotDf = data.frame()
    markPeakFiles = Sys.glob(file.path(baseDir, sprintf("roadmap/peaks/narrowPeak/*%s*.narrowPeak.gz", peakType)))
    for (markPeakFile in markPeakFiles) {
        print(markPeakFile)
        eid = gsub("(E\\d+)-.*", "\\1", basename(markPeakFile), perl = T)
        markPeakDf = read.delim(gzfile(markPeakFile), header = F, as.is = T)[, c(1:3)]
        colnames(markPeakDf) = c("space", "start", "end")
        markPeakDf = markPeakDf[markPeakDf$space %in% cchrs,]
        markPeakDf$eid = eid
        markPeakTotDf = rbind(markPeakTotDf, markPeakDf)
    }
    markPeakTotDf$space = factor(markPeakTotDf$space, levels = cchrs)
    markPeakTotDf = markPeakTotDf[order(markPeakTotDf$space, markPeakTotDf$start, markPeakTotDf$end),]
    markPeakTotGr = GRanges(seqnames = Rle(markPeakTotDf$space),
                            IRanges(start = markPeakTotDf$start,
                                    end = markPeakTotDf$end),
                            strand = "*",
                            eid = markPeakTotDf$eid)

    markPeakTotReducedGr = reduce(markPeakTotGr)
    markPeakTotReducedGr$ecnt = countOverlaps(markPeakTotReducedGr, markPeakTotGr)
    #ol = findOverlaps(markPeakTotReducedGr, markPeakTotGr)
    #markPeakTotReducedGr$eids = unlist(mclapply(unique(queryHits(ol)),
                                        #function (i) {
                                            #paste(markPeakTotGr[subjectHits(ol[queryHits(ol) == i])]$eid, collapse = ",")
                                        #}, mc.cores = 10))

    cPeakType = capitalize(peakType)
    #markPeakTotReducedDf = as.data.frame(markPeakTotReducedGr)[, c(1,2,3,6,7)]
    markPeakTotReducedDf = as.data.frame(markPeakTotReducedGr)[, c(1,2,3,6)]
    markPeakTotReducedFile = file.path(baseDir, sprintf("roadmap/peaks/roadmap%sNarrowPeaks.bed", cPeakType))
    write.table(markPeakTotReducedDf, markPeakTotReducedFile, sep = "\t", row.names=F, col.names=F, quote=FALSE)
}


