rm(list=ls())

##
## Load libraries
##
require(doMC)
registerDoMC(8)

require(proxy)
require(sqldf)
require(gdata)
require(gtools)
require(annotate)
require(snpStats)
require(data.table)
require(Biostrings)
require(GenomicRanges)
require(org.Hs.eg.db)
require(BSgenome.Hsapiens.UCSC.hg19)


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
rootDir = "/home/sl279"
baseDir = file.path(rootDir, "BiO/Research/NoncoDiver")
vcfDir = file.path(baseDir, "vcf")
panVcfDir = file.path(vcfDir, "pancan")
scntFiles = mixedsort(Sys.glob(file.path(panVcfDir, "*.scnt.txt")))

chrs = c(1:22, "X", "Y")
chrs = factor(chrs, level=chrs)
cchrs = sapply(chrs, function(x) paste("chr", x, sep=""))
cchrs = factor(cchrs, level=cchrs)
distCut = 100
rscntPerSloci = 2

foreach(scntFile=scntFiles) %dopar% {
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

    #length(scntGrp)

    scntGrpStartPos = lapply(scntGrp, function(x) { min(scntDf[x, "pos"]) })
    scntGrpStartPos = unlist(scntGrpStartPos)

    #length(scntGrpStartPos)

    scntGrpEndPos = lapply(scntGrp, function(x) { max(scntDf[x, "pos"]) })
    scntGrpEndPos = unlist(scntGrpEndPos)

    scntGrpSnvLoci = sapply(1:length(scntGrp), function(i) { nrow(scntDf[scntGrp[[i]],]) })
    scntGrpSnvLoci = unlist(scntGrpSnvLoci)

    scntGrpSnvCnt = sapply(1:length(scntGrp), function(i) { sum( as.numeric(unlist(sapply(scntDf[scntGrp[[i]], "scnt"], function(x) strsplit(x, ",", fixed = T)[[1]])))) } )
    scntGrpSnvCnt = unlist(scntGrpSnvCnt)

    #length(scntGrpSnvCnt)
    #summary(scntGrpSnvCnt)

    scntGrpSnvRate = sapply(1:length(scntGrp), function(x) { scntGrpSnvCnt[x] / (max(scntDf[scntGrp[[x]], "pos"]) - min(scntDf[scntGrp[[x]], "pos"]) + 1) })
    scntGrpSnvRate = unlist(scntGrpSnvRate)

    #length(scntGrpSnvRate)
    #summary(scntGrpSnvRate)

    scntGrpStartSids = lapply(scntGrp, function(x) { unique(as.vector(unlist(sapply(as.vector(unlist(sapply(scntDf[x, "sids"],
                                                                          function(x) strsplit(x, ",", fixed = T)[[1]]))),
                                                                   function(y) strsplit(y, ";", fixed = T)[[1]])))) })
    scntGrpStartSids = lapply(scntGrpStartSids, function(x) { paste(x, collapse = ",") } )
    scntGrpStartSids = unlist(scntGrpStartSids)

    #head(scntGrpStartSids)
    #length(scntGrpStartSids)
    #summary(elementLengths(scntGrpStartSids))

    scntGrpNrSnvCnt = as.vector(sapply(scntGrpStartSids, function(x) length(strsplit(x, ",", fixed = T)[[1]])))
    #summary(scntGrpNrSnvCnt)

    scntGrpSeqs = getSeq(BSgenome.Hsapiens.UCSC.hg19,
                         GRanges(seqnames = Rle(gsub("(.*)", "chr\\1", rep(scntDf$chrom[1], length(scntGrp)), perl = T)),
                                 ranges = IRanges(start = scntGrpStartPos, end = scntGrpEndPos),
                                 strand = "+"),
                         as.character = T)

    scntGrpNrSnvRef = sapply(1:length(scntGrp), function(x) {
                                 if (length(scntGrp[[x]]) == 1) {
                                     return(scntDf[scntGrp[[x]], "ref"])
                                 } else {
                                     return(scntGrpSeqs[x])
                                 }
                         })

    scntGrpNrSnvAlt = sapply(1:length(scntGrp), function(x) {
                                 if (length(scntGrp[[x]]) == 1) {
                                     return(scntDf[scntGrp[[x]], "alt"])
                                 } else {
                                     return("-")
                                 }
                         })

    scntGrpCaddRaw = sapply(1:length(scntGrp), function(x) { mean(as.numeric(sapply(scntDf[x, "cadd_raw"], function(x) strsplit(x, ",", fixed = T)[[1]]))) } )
    scntGrpCaddPhred = sapply(1:length(scntGrp), function(x) { mean(as.numeric(sapply(scntDf[x, "cadd_phred"], function(x) strsplit(x, ",", fixed = T)[[1]]))) } )

    scntGrpDf = data.frame(space = rep(scntDf[1, "chrom"], length(scntGrp)),
                           start = scntGrpStartPos,
                           end = scntGrpEndPos,
                           width = scntGrpEndPos - scntGrpStartPos + 1,
                           ref = scntGrpNrSnvRef,
                           alt = scntGrpNrSnvAlt,
                           allele = paste(scntGrpNrSnvRef, scntGrpNrSnvAlt, sep = "/"),
                           strand = "+",
                           sloci = scntGrpSnvLoci,
                           rscnt = scntGrpSnvCnt,
                           nrscnt = scntGrpNrSnvCnt,
                           rscnt_rate = scntGrpSnvRate,
                           rscnt_per_sloci = scntGrpSnvCnt / scntGrpSnvLoci,
                           sids = scntGrpStartSids,
                           avg_cadd_raw = scntGrpCaddRaw,
                           avg_cadd_phred = scntGrpCaddPhred)

    scntGrpFile = gsub("scnt", sprintf("hotspot%d", distCut), scntFile)
    write.table(scntGrpDf, scntGrpFile, row.names = F, col.names = T, sep = "\t", quote = F)

    scntGrpVepDf = scntGrpDf
    scntGrpVepFile = gsub(sprintf("hotspot%d", distCut), sprintf("hotspot%d.vep_in", distCut), scntGrpFile)
    write.table(scntGrpVepDf[, c("space", "start", "end", "allele", "strand")],
                scntGrpVepFile, row.names = F, col.names = F, sep = "\t", quote = F)

    scntGrpVepSigDf = subset(scntGrpVepDf, rscnt_per_sloci > 1)
    scntGrpVepSigFile = gsub(sprintf("hotspot%d", distCut), sprintf("hotspot%d.sig%d", distCut, rscntPerSloci), scntGrpVepFile)
    write.table(scntGrpVepSigDf[, c("space", "start", "end", "allele", "strand")],
                scntGrpVepSigFile, row.names = F, col.names = F, sep = "\t", quote = F)
}

##
hotspotChrVepSigFiles = mixedsort(Sys.glob(file.path(panVcfDir, sprintf("*hotspot%d.sig%d.vep_in.txt", distCut, rscntPerSloci))))
#hotspotChrVepSigFiles = mixedsort(Sys.glob(file.path(panVcfDir, sprintf("*hotspot%d.vep_in.txt", distCut))))

for (hotspotChrVepSigFile in hotspotChrVepSigFiles) {
    print(hotspotChrVepSigFile)
    hotspotChrVepSigOutFile = gsub("vep_in", "vep_out", hotspotChrVepSigFile)
    hotspotChrVepSigLsfOutFile = gsub("txt", "lsfout", hotspotChrVepSigOutFile)
    system(sprintf("
bsub -g /nd/hotspot/vep \\
    -q i2b2_12h -W 12:0 \\
    -n 8 -R \"span[hosts=1]\" \\
    -o %s \\
    perl /home/sl279/vep/variant_effect_predictor.pl \\
        --fork 8 \\
        --force_overwrite \\
        --offline \\
        --no_stats \\
        --everything \\
        --check_existing \\
        --total_length \\
        --allele_number \\
        --no_escape \\
        --gencode_basic \\
        --assembly GRCh37 \\
        --dir /home/sl279/.vep \\
        --plugin UpDownDistance,1000000,1000000 \\
        --fasta /home/sl279/.vep/homo_sapiens/76_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/database/gap.bed.gz,gap,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/database/genomicSuperDups.bed.gz,genomicSuperDups,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/database/gwasCatalog.bed.gz,gwasCatalog,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/database/targetScanS.bed.gz,targetScanS,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/database/tfbsConsSites.bed.gz,tfbsConsSites,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/database/wgRna.bed.gz,wgRna,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/database/phastConsElements46way.bed.gz,phastConsElements46way,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/database/phastConsElements100way.bed.gz,phastConsElements100way,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/database/rmsk.bed.gz,rmsk,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/database/nestedRepeats.bed.gz,nestedRepeats,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/database/simpleRepeat.bed.gz,simpleRepeat,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/database/microsat.bed.gz,microsat,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/encode/wgEncodeAwgDnaseMasterSites.bed.gz,wgEncodeAwgDnaseMasterSites,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/encode/wgEncodeRegDnaseClusteredV3.bed.gz,wgEncodeRegDnaseClusteredV3,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/fantom5/fantom5PermissiveEnhancers.bed.gz,fantom5PermissiveEnhancers,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/dbsuper/all/allDbSuperEnhancers.bed.gz,allDbSuperEnhancers,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/encode/wgEncodeRegTfbsClusteredWithCellsV3.bed.gz,wgEncodeRegTfbsClusteredWithCellsV3,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/encode/wgEncodeAwgSegmentationCombinedGm12878.bed.gz,wgEncodeAwgSegmentationCombinedGm12878,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/encode/wgEncodeAwgSegmentationCombinedH1hesc.bed.gz,wgEncodeAwgSegmentationChromhmmH1hesc,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/encode/wgEncodeAwgSegmentationCombinedHelas3.bed.gz,wgEncodeAwgSegmentationCombinedHelas3,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/encode/wgEncodeAwgSegmentationCombinedHepg2.bed.gz,wgEncodeAwgSegmentationCombinedHepg2,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/encode/wgEncodeAwgSegmentationCombinedHuvec.bed.gz,wgEncodeAwgSegmentationCombinedHuvec,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/encode/wgEncodeAwgSegmentationCombinedK562.bed.gz,wgEncodeAwgSegmentationCombinedK562,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/encode/distalDhsToPromoterDhs.bed.gz,distalDhsToPromoterDhs,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/encode/dhsToGeneExpression.bed.gz,dhsToGeneExpression,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/fantom5/fantom5EnhancerTssAssociations.bed.gz,fantom5EnhancerTssAssociations,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/dbsuper/all/allDbSuperEnhancerGeneAssociations.bed.gz,allDbSuperEnhancerGeneAssociations,bed,overlap,0 \\
        --input_file %s \\
        --output_file %s
", hotspotChrVepSigLsfOutFile, hotspotChrVepSigFile, hotspotChrVepSigOutFile), intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE, wait = F)
}


##
hotspotChrFiles = mixedsort(Sys.glob(file.path(panVcfDir, sprintf("TCGA_16_Cancer_Types.wgs.somatic.chr*.sanitized.hotspot%d.txt", distCut))))
hotspotDf = data.frame()
for (hotspotChrFile in hotspotChrFiles) {
    cat(sprintf("Reading %s ...\n", hotspotChrFile))
    hotspotChrDf = read.delim(hotspotChrFile, header = T, as.is = T)
    hotspotChrDf$rscnt_per_sloci = hotspotChrDf$rscnt / hotspotChrDf$sloci
    hotspotDf = rbind(hotspotDf, hotspotChrDf)
    hotspotChrSigDf = subset(hotspotChrDf, rscnt_per_sloci > 1)
    hotspotChrSigFile = gsub(sprintf("hotspot%d", distCut), sprintf("hotspot%d.sig%d.vep_in", distCut, rscntPerSloci), hotspotChrFile)
    #cat(sprintf("Writing %s ...\n", hotspotChrSigFile))
    #write.table(hotspotChrSigDf[, c("space", "start", "end", "allele", "strand")],
                #hotspotChrSigFile, row.names = F, col.names = F, sep = "\t", quote = F)
}

hotspotFile = file.path(panVcfDir, sprintf("TCGA_16_Cancer_Types.wgs.somatic.sanitized.hotspot%d.txt", distCut))
write.table(hotspotDf, hotspotFile, row.names = F, col.names = T, sep = "\t", quote = F)

names(hotspotDf)
head(hotspotDf)
class(hotspotDf$rscnt)
class(hotspotDf$sloci)
summary(hotspotDf$rscnt)
summary(hotspotDf$nrscnt)
summary(hotspotDf$rscnt_per_sloci)
summary(hotspotDf$avg_cadd_raw)
summary(hotspotDf$avg_cadd_phred)

hotspotDf = hotspotDf[order(-hotspotDf$rscnt_per_sloci),]
hotspotSigDf = subset(hotspotDf, rscnt_per_sloci > 1)

##
## Expand vep_out files using Ruby script
##


##
## Annotate significant hotspots with VEP output
##
hotspotSigVepOutChrFiles = mixedsort(Sys.glob(file.path(panVcfDir, sprintf("TCGA_16_Cancer_Types.wgs.somatic.chr*.hotspot%d.sig%d.vep_out.exp.txt", distCut, rscntPerSloci))))
hotspotSigVepOutDf = data.frame()
for (hotspotSigVepOutChrFile in hotspotSigVepOutChrFiles) {
    cat(sprintf("Reading %s ...\n", hotspotSigVepOutChrFile))
    hotspotSigVepOutChrDf = read.delim(hotspotSigVepOutChrFile, header = T, as.is = T)
    hotspotSigVepOutChrDf = subset(hotspotSigVepOutChrDf, genomicSuperDups == "" & !grepl("trf", simpleRepeat))
    hotspotSigVepOutDf = rbind(hotspotSigVepOutDf, hotspotSigVepOutChrDf)
}

hotspotSigVepOutFile = file.path(panVcfDir, sprintf("TCGA_16_Cancer_Types.wgs.somatic.sanitized.hotspot%d.sig%d.vep_out.exp.txt", distCut, rscntPerSloci))
write.table(hotspotSigVepOutDf, hotspotSigVepOutFile, row.names = F, col.names = T, sep = "\t", quote = F)
#hotspotSigVepOutDf = read.delim(hotspotSigVepOutFile, header = T, as.is = T)
head(hotspotSigVepOutDf)
tail(hotspotSigVepOutDf)

hotspotSigVepOutDf[, "space"] = gsub("(^\\S+?):.*", "\\1", hotspotSigVepOutDf$Location, perl = T)
hotspotSigVepOutDf[, "start"] = sapply(hotspotSigVepOutDf$Location, function(x) { strsplit(strsplit(x, ":", fixed = T)[[1]][2], "-", fixed = T)[[1]][1] })
hotspotSigVepOutDf[, "end"] = sapply(hotspotSigVepOutDf$Location, function(x) { strsplit(strsplit(x, ":", fixed = T)[[1]][2], "-", fixed = T)[[1]][2] })
hotspotSigVepOutDf[is.na(hotspotSigVepOutDf$end), "end"] = hotspotSigVepOutDf[is.na(hotspotSigVepOutDf$end), "start"]
hotspotSigVepOutMergedDf = merge(hotspotSigVepOutDf, hotspotSigDf, by = c("space", "start", "end"), all.x = TRUE)
hotspotSigVepOutMergedDf = hotspotSigVepOutMergedDf[order(-hotspotSigVepOutMergedDf$rscnt_per_sloci),]
head(hotspotSigVepOutMergedDf, 20)

hotspotSigVepOutFile2 = file.path(panVcfDir, sprintf("TCGA_16_Cancer_Types.wgs.somatic.sanitized.hotspot%d.sig%d.vep_out.exp.with_scnt.txt", distCut, rscntPerSloci))
write.table(hotspotSigVepOutMergedDf, hotspotSigVepOutFile2, row.names = F, col.names = T, sep = "\t", quote = F)

save.image()
