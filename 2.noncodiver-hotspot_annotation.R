rm(list=ls())

##
## Load libraries
##
require(doMC)
registerDoMC(12)

require(hash)
require(Rcpp)
require(RcppArmadillo)
require(proxy)
require(sqldf)
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
panVcfDir = file.path(vcfDir, "pancan")
gdacDir = file.path(baseDir, "gdac")
gdacStdDataDir = file.path(gdacDir, "stddata__2015_02_04")
gdacAnlDataDir = file.path(gdacDir, "analyses__2014_10_17")
figDir = file.path(baseDir, "figure")

chrs = c(1:22, "X", "Y")
chrs = factor(chrs, level=chrs)
cchrs = sapply(chrs, function(x) paste("chr", x, sep=""))
cchrs = factor(cchrs, level=cchrs)
distCut = 100
margin = 50
permNo = 100
rscntPerSloci = 2
baseFontSize = 15

# Load chrom size information
hg19File = file.path(baseDir, "ucsc/database/hg19.genome")
hg19Df = read.delim(hg19File, header=T, as.is=T)
hg19Df = hg19Df[hg19Df$chrom %in% chrs,]
hg19Df$chrom = factor(hg19Df$chrom, levels=chrs)
hg19Df = hg19Df[order(hg19Df$chrom),]

# Cluster mutations
scntFiles = mixedsort(Sys.glob(file.path(panVcfDir, "*.scnt.txt")))

#foreach(scntFile=scntFiles) %dopar% {
for (scntFile in scntFiles) {
    cat(sprintf("Processing %s ...\n", scntFile))
    scntDf = read.delim(scntFile, header = F, as.is = T)
    colnames(scntDf) = c("chrom", "pos", "ref", "alt", "scnt", "sids", "cadd_raw", "cadd_phred")
    #scntDf$tscnt = sapply(scntDf$scnt, function(x) { sum( as.numeric(unlist(sapply(x, function(x) strsplit(x, ",", fixed = T)[[1]])))) } )
    distPos <- foreach(i=2:nrow(scntDf), .combine = c) %dopar% {
        return(scntDf[i,]$pos - scntDf[i-1,]$pos)
    }
    summary(distPos)
    hist(distPos)

    # Calculate total number of SNVs for each patient
    allSids = unique(unlist(sapply(scntDf$sids, function(x) { unlist(strsplit(strsplit(x, ",")[[1]], ";")) })))
    chrScnt = list()
    for (sid in allSids) {
        chrScnt[[sid]] = length(grep(sid, scntDf$sids))
    }

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

    grpDistPos <- foreach(i=2:length(scntGrp), .combine = c) %dopar% {
        prevGrp = scntGrp[[i-1]]
        currGrp = scntGrp[[i]]
        return(min(scntDf[currGrp,]$pos) - max(scntDf[prevGrp,]$pos))
    }
    summary(grpDistPos)

    #scntGrpSids = lapply(scntGrp, function(x) { unique(as.vector(unlist(sapply(as.vector(unlist(sapply(scntDf[x, "sids"],
                                                                                                       #function(x) strsplit(x, ",", fixed = T)[[1]]))),
                                                                               #function(y) strsplit(y, ";", fixed = T)[[1]])))) })

    scntGrpSids = lapply(scntGrp, function(x) { unique(unlist(sapply(scntDf[x, "sids"],
                                                                     function(x) { unlist(strsplit(strsplit(x, ",")[[1]], ";")) }))) })

    #scntGrpSids = sort(unique(unlist(scntGrpSids)))
    #scntGrpSids2 = sort(unique(unlist(scntGrpSids2)))
    #length(scntGrpSids)
    #length(scntGrpSids2)
    #identical(scntGrpSids, scntGrpSids2)

    scntGrpSids = lapply(scntGrpSids, function(x) { paste(x, collapse = ",") } )
    scntGrpSids = unlist(scntGrpSids)

    scntGrpSidCnt = as.vector(sapply(scntGrpSids, function(x) length(strsplit(x, ",", fixed = T)[[1]])))
    summary(scntGrpSidCnt)

    scntSelGrp = scntGrp[scntGrpSidCnt > 1]
    length(scntSelGrp)

    scntSelGrpSids = lapply(scntSelGrp, function(x) { unique(as.vector(unlist(sapply(as.vector(unlist(sapply(scntDf[x, "sids"],
                                                                                                       function(x) strsplit(x, ",", fixed = T)[[1]]))),
                                                                               function(y) strsplit(y, ";", fixed = T)[[1]])))) })
    scntSelGrpSids = lapply(scntSelGrpSids, function(x) { paste(x, collapse = ",") } )
    scntSelGrpSids = unlist(scntSelGrpSids)

    scntSelGrpSidCnt = as.vector(sapply(scntSelGrpSids, function(x) length(strsplit(x, ",", fixed = T)[[1]])))

    scntSelGrpStartPos = lapply(scntSelGrp, function(x) { min(scntDf[x, "pos"]) - margin })
    scntSelGrpStartPos = unlist(scntSelGrpStartPos)

    scntSelGrpEndPos = lapply(scntSelGrp, function(x) { max(scntDf[x, "pos"]) + margin })
    scntSelGrpEndPos = unlist(scntSelGrpEndPos)

    scntSelGrpSnvLoci = sapply(1:length(scntSelGrp), function(i) { nrow(scntDf[scntSelGrp[[i]],]) })
    scntSelGrpSnvLoci = unlist(scntSelGrpSnvLoci)

    scntSelGrpSnvCnt = sapply(1:length(scntSelGrp), function(i) { sum( as.numeric(unlist(sapply(scntDf[scntSelGrp[[i]], "scnt"], function(x) strsplit(x, ",", fixed = T)[[1]])))) } )
    scntSelGrpSnvCnt = unlist(scntSelGrpSnvCnt)

    scntSelGrpSnvRate = sapply(1:length(scntSelGrp), function(x) { scntSelGrpSnvCnt[x] / (max(scntDf[scntSelGrp[[x]], "pos"]) - min(scntDf[scntSelGrp[[x]], "pos"]) + 1) })
    scntSelGrpSnvRate = unlist(scntSelGrpSnvRate)

    chrL = hg19Df[hg19Df$chrom == scntDf$chrom[1],]$size

    #cppFunction('int sumSampleMutCnt(int grpB, int grpE, int chrL, std::vector< std::string > allSids, List chrScnt, Function f) {
                #std::vector<int> sMcnts;
                #for (std::vector< std::string >::iterator it = allSids.begin(); it != allSids.end(); ++it) {
                    #IntegerVector mutations = f(chrL, chrScnt[*it]);
                    #int sMcnt = 0;
                    #for (IntegerVector::iterator it2 = mutations.begin(); it2 != mutations.end(); ++it2) {
                        #if (*it2 >= grpB && *it2 <= grpE) {
                            #sMcnt += 1;
                        #}
                    #}
                    #sMcnts.push_back(sMcnt);
                #}
                #return(std::accumulate(sMcnts.begin(), sMcnts.end(), 0)); }')

    scntSelGrpPvalues = mclapply(1:length(scntSelGrp),
                               function(i) {
                                   grpB = scntSelGrpStartPos[i]
                                   grpE = scntSelGrpEndPos[i]
                                   grpL = grpE - grpB + 1
                                   grpS = scntSelGrpSidCnt[i]
                                   smpB = grpB - 5e3
                                   smpE = grpE + 5e3
                                   smpB = ifelse(smpB < 0, 0, smpB)
                                   smpE = ifelse(smpE > chrL, chrL, smpE)
                                   smpL = smpE - smpB + 1
                                   smpS = length(unique(unlist(sapply(unlist(sapply(subset(scntDf, pos >= smpB & pos <= smpE)$sids,
                                                                                    function(x) { strsplit(x, ",")[[1]] })),
                                                                             function(y) { strsplit(y, ";")[[1]] }))))
                                   lambda = smpS / smpL
                                   grpPoiP = poisson.test(grpS, lambda * grpL, alternative = "greater")$p.value

                                   #rndS <- foreach(i=1:100000, .combine = c) %dopar% {
                                       #sum(sapply(allSids,
                                                  #function(y) {
                                                      ##smpS = length(grep(y, subset(scntDf, pos >= smpB & pos <= smpE)$sids))
                                                      ##mutP = sample(c(smpB:smpE), smpS)
                                                      #mutP = sample.int(chrL, chrScnt[[y]]);
                                                      #return(length(mutP[mutP >= grpB & mutP <= grpE]) > 0)
                                                  #}))
                                   #}
                                   #table(rndS)
                                   #hist(rndS)
                                   #plot(density(rndS))

                                   #system.time(
                                   #rndS <- foreach(i=1:permNo, .combine = c) %dopar% {
                                       #sumSampleMutCnt(grpB, grpE, chrL, allSids, chrScnt, sample.int)
                                   #}
                                   #)

                                   #rndS <- calcPermMutCnt(permNo, grpB, grpE, chrL, allSids, chrScnt, sample.int)
                                   #grpPermP = (sum(rndS >= grpS) + 1) / (permNo + 1)
                                   #grpHp = prod(sapply(allSids, function(y) dhyper(length(grep(y, scntDf[x, "sids"])), grpL, chrL - grpL, chrScnt[[y]])))
                                   #grpFp = Fisher.test(sapply(allSids, function(y) dhyper(length(grep(y, scntDf[x, "sids"])), grpL, chrL - grpL, chrScnt[[y]])))[2]
                                   #grpBp = Stouffer.test(sapply(allSids, function(y) dhyper(length(grep(y, scntDf[x, "sids"])), grpL, chrL - grpL, chrScnt[[y]])))[2]
                                   #cat(sprintf("(%d/%d) Hotspot %s:%d-%d, permutation p-value = %f, multiplied hypergeometric p-value = %f, Fisher's combined hypergeometric p-value = %f, Stouffer's combined p-value = %f\n", 
                                               #i, length(scntSelGrp), scntDf$chrom[1], grpB, grpE, grpPp, grpHp, grpFp, grpBp))
                                   cat(sprintf("(%d/%d) Hotspot chr%s:%d-%d, Observed sample count: %d, Poisson p-value = %f\n", i, length(scntSelGrp), scntDf$chrom[1], grpB, grpE, grpS, grpPoiP))
                                   return(grpPoiP)
                               }, mc.cores = 12)
                             
    scntSelGrpPvalues = unlist(scntSelGrpPvalues)
    #summary(scntSelGrpPvalues)
    #hist(scntSelGrpPvalues)

    scntSelGrpFdrs = p.adjust(scntSelGrpPvalues, method = "fdr")
    #summary(scntSelGrpFdrs)
    #hist(scntSelGrpFdrs)

    scntSelGrpSeqs = getSeq(BSgenome.Hsapiens.UCSC.hg19,
                         GRanges(seqnames = Rle(gsub("(.*)", "chr\\1", rep(scntDf$chrom[1], length(scntSelGrp)), perl = T)),
                                 ranges = IRanges(start = scntSelGrpStartPos, end = scntSelGrpEndPos),
                                 strand = "+"),
                         as.character = T)

    scntSelGrpNrSidRef = sapply(1:length(scntSelGrp), function(x) {
                                 if (length(scntSelGrp[[x]]) == 1) {
                                     return(scntDf[scntSelGrp[[x]], "ref"])
                                 } else {
                                     return(scntSelGrpSeqs[x])
                                 }
                         })

    scntSelGrpNrSidAlt = sapply(1:length(scntSelGrp), function(x) {
                                 if (length(scntSelGrp[[x]]) == 1) {
                                     return(scntDf[scntSelGrp[[x]], "alt"])
                                 } else {
                                     return("-")
                                 }
                         })

    scntSelGrpCaddRaw = sapply(1:length(scntSelGrp), function(x) { mean(as.numeric(sapply(scntDf[x, "cadd_raw"], function(x) strsplit(x, ",", fixed = T)[[1]]))) } )
    scntSelGrpCaddPhred = sapply(1:length(scntSelGrp), function(x) { mean(as.numeric(sapply(scntDf[x, "cadd_phred"], function(x) strsplit(x, ",", fixed = T)[[1]]))) } )

    scntSelGrpDf = data.frame(space = rep(scntDf[1, "chrom"], length(scntSelGrp)),
                           start = scntSelGrpStartPos,
                           end = scntSelGrpEndPos,
                           width = scntSelGrpEndPos - scntSelGrpStartPos + 1,
                           ref = scntSelGrpNrSidRef,
                           alt = scntSelGrpNrSidAlt,
                           allele = paste(scntSelGrpNrSidRef, scntSelGrpNrSidAlt, sep = "/"),
                           strand = "+",
                           sloci = scntSelGrpSnvLoci,
                           rscnt = scntSelGrpSnvCnt,
                           nrscnt = scntSelGrpSidCnt,
                           rscnt_rate = scntSelGrpSnvRate,
                           rscnt_per_sloci = scntSelGrpSnvCnt / scntSelGrpSnvLoci,
                           p_value = scntSelGrpPvalues,
                           adj_p_value = scntSelGrpFdrs,
                           sids = scntSelGrpSids,
                           avg_cadd_raw = scntSelGrpCaddRaw,
                           avg_cadd_phred = scntSelGrpCaddPhred)

    scntSelGrpDf.2 = scntSelGrpDf[order(scntSelGrpDf$adj_p_value),]
    head(scntSelGrpDf.2)

    scntSelGrpFile = gsub("scnt", sprintf("hotspot%d", distCut), scntFile)
    write.table(scntSelGrpDf, scntSelGrpFile, row.names = F, col.names = T, sep = "\t", quote = F)

    scntSelGrpVepDf = scntSelGrpDf
    scntSelGrpVepFile = gsub(sprintf("hotspot%d", distCut), sprintf("hotspot%d.vep_in", distCut), scntSelGrpFile)
    write.table(scntSelGrpVepDf[, c("space", "start", "end", "allele", "strand")],
                scntSelGrpVepFile, row.names = F, col.names = F, sep = "\t", quote = F)

}

##
hotspotChrVepSigFiles = mixedsort(Sys.glob(file.path(panVcfDir, sprintf("*hotspot%d.vep_in.txt", distCut))))

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
}

hotspotFile = file.path(panVcfDir, sprintf("TCGA_16_Cancer_Types.wgs.somatic.sanitized.hotspot%d.txt", distCut))
write.table(hotspotDf, hotspotFile, row.names = F, col.names = T, sep = "\t", quote = F)

#names(hotspotDf)
#head(hotspotDf)
#class(hotspotDf$rscnt)
#class(hotspotDf$sloci)
#summary(hotspotDf$rscnt)
#summary(hotspotDf$nrscnt)
#summary(hotspotDf$rscnt_per_sloci)
#summary(hotspotDf$avg_cadd_raw)
#summary(hotspotDf$avg_cadd_phred)

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

#hotspotSigVepOutFile = file.path(panVcfDir, sprintf("TCGA_16_Cancer_Types.wgs.somatic.sanitized.hotspot%d.sig%d.vep_out.exp.txt", distCut, rscntPerSloci))
#write.table(hotspotSigVepOutDf, hotspotSigVepOutFile, row.names = F, col.names = T, sep = "\t", quote = F)
#hotspotSigVepOutDf = read.delim(hotspotSigVepOutFile, header = T, as.is = T)
#head(hotspotSigVepOutDf)
#tail(hotspotSigVepOutDf)

hotspotSigVepOutDf[, "space"] = gsub("(^\\S+?):.*", "\\1", hotspotSigVepOutDf$Location, perl = T)
hotspotSigVepOutDf[, "start"] = sapply(hotspotSigVepOutDf$Location, function(x) { strsplit(strsplit(x, ":", fixed = T)[[1]][2], "-", fixed = T)[[1]][1] })
hotspotSigVepOutDf[, "end"] = sapply(hotspotSigVepOutDf$Location, function(x) { strsplit(strsplit(x, ":", fixed = T)[[1]][2], "-", fixed = T)[[1]][2] })
hotspotSigVepOutDf[is.na(hotspotSigVepOutDf$end), "end"] = hotspotSigVepOutDf[is.na(hotspotSigVepOutDf$end), "start"]

hotspotSigVepOutMergedDf = merge(hotspotSigVepOutDf, hotspotSigDf, by = c("space", "start", "end"), all.x = TRUE)
hotspotSigVepOutMergedDf = hotspotSigVepOutMergedDf[order(-hotspotSigVepOutMergedDf$nrscnt),]
head(hotspotSigVepOutMergedDf)

hotspotSigVepOutFile2 = file.path(panVcfDir, sprintf("TCGA_16_Cancer_Types.wgs.somatic.sanitized.hotspot%d.sig%d.vep_out.exp.with_scnt.txt", distCut, rscntPerSloci))
write.table(hotspotSigVepOutMergedDf, hotspotSigVepOutFile2, row.names = F, col.names = T, sep = "\t", quote = F)

##
## Test mRNA expression changes 
##

### Load RefSeq gene info table (ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz)
#geneInfoFile = file.path(baseDir, "entrez/gene_info_human.txtd")
#geneInfoDt = fread(geneInfoFile, header = F)
#geneInfoDt = geneInfoDt[, c(2:5), with = F]
#setnames(geneInfoDt, c("GeneID", "Symbol", "LocusTag", "Synonyms"))
#geneInfoExtDt = geneInfoDt[, list(Synonym = unlist(strsplit(Synonyms, "|", fixed = T), use.names = FALSE)), by = "GeneID,Symbol,LocusTag"]
#geneInfoExtDf = as.data.frame(geneInfoExtDt)

## Load mRNA expression data from GDAC
rsemFiles = Sys.glob(file.path(gdacStdDataDir, "*/*/*/*RSEM_all.txt"))
rsemPanDf = data.frame()
rsemLog2FoldPanDf = data.frame()
for (i in 1:length(rsemFiles)) {
    rsemFile = rsemFiles[i]
    cat(sprintf("Processing %s ...\n", rsemFile))
    rsemDf = read.delim(rsemFile, header = T, as.is = T)
    colnames(rsemDf)[1] = "id"
    colnames(rsemDf)[2:ncol(rsemDf)] = as.vector(sapply(colnames(rsemDf)[2:ncol(rsemDf)], function(x) paste(strsplit(x, ".", fixed = T)[[1]][2:4], collapse = "__")))
    rsemDf$symbol = sapply(rsemDf$id, function(x) strsplit(x, "|", fixed = T)[[1]][1])
    rsemDf$entrez = sapply(rsemDf$id, function(x) strsplit(x, "|", fixed = T)[[1]][2])

    rsemLog2Df = log2(rsemDf[, grep("__", colnames(rsemDf))] + 1)
    rsemLog2FoldDf = rsemLog2Df - rowMeans(rsemLog2Df)
    rsemLog2FoldDf = data.frame(id = rsemDf$id, symbol = rsemDf$symbol, entrez = rsemDf$entrez, rsemLog2FoldDf)
    if (i == 1) {
        rsemPanDf = rsemDf
        rsemLog2FoldPanDf = rsemLog2FoldDf
    } else {
        rsemPanDf = cbind(rsemPanDf, rsemDf[, grep("__", colnames(rsemDf))])
        rsemLog2FoldPanDf = cbind(rsemLog2FoldPanDf, rsemLog2FoldDf[, grep("__", colnames(rsemLog2FoldDf))])
    }
}

### Update gene names using Entrez gene ID
#newGeneNames = as.vector(unlist(mclapply(rsemPanDf$id, function(x) {
    #rSymbol = strsplit(x, "|", fixed = T)[[1]][1]
    #rEntrez = strsplit(x, "|", fixed = T)[[1]][2]
    #rEntrez2Symbol = getSYMBOL(rEntrez, data='org.Hs.eg')
    #nSymbol = ifelse(is.na(rEntrez2Symbol), rSymbol, rEntrez2Symbol)
    #return(nSymbol)
#}, mc.cores = 8)))


##
## Load SCNA data for filter
##
cnvByGeneFiles = Sys.glob(file.path(gdacAnlDataDir, "*/*/*/all_data_by_genes.txt"))
cnvByGeneDf = data.frame()
for (i in 1:length(cnvByGeneFiles)) {
    cnvByGeneFile = cnvByGeneFiles[i]
    cat(sprintf("Reading %s ...\n", cnvByGeneFile))
    cnvByGeneTmpDf = read.delim(cnvByGeneFile, header = T, as.is = T)
    cnvByGeneTmpDf = cnvByGeneTmpDf[, c(-2,-3)]
    colnames(cnvByGeneTmpDf)[2:ncol(cnvByGeneTmpDf)] = as.vector(sapply(colnames(cnvByGeneTmpDf)[2:ncol(cnvByGeneTmpDf)],
                                                                        function(x) {
                                                                            elems = strsplit(x, ".", fixed = T)[[1]][2:4]
                                                                            elems[3] = substring(elems[3], 1, 2)
                                                                            paste(elems, collapse = "__")
                                                                        }))
    if (i == 1) {
        cnvByGeneDf = cnvByGeneTmpDf
    } else {
        cnvByGeneDf = cbind(cnvByGeneDf, cnvByGeneTmpDf[, grep("__", colnames(cnvByGeneTmpDf))])
    }
}
rownames(cnvByGeneDf) = cnvByGeneDf$Gene.Symbol
cnvByGeneDf$Gene.Symbol = NULL

##
#hotspotSigVepOutMergedDf$allele = NULL
#hotspotSigVepOutMergedDf$strand = NULL
hotspotSigGenesDf = sqldf('SELECT space, start, end, SYMBOL, sids
                          FROM hotspotSigVepOutMergedDf
                          WHERE nrscnt > 10
                          GROUP BY space, start, end, SYMBOL
                          ORDER BY nrscnt DESC')

#head(hotspotSigGenesDf)
#tail(hotspotSigGenesDf)
#nrow(hotspotSigGenesDf)

#for (i in 1:nrow(hotspotSigGenesDf)) {
hotspotSigGenesNewDf <- foreach (i=1:nrow(hotspotSigGenesDf), .combine=rbind) %dopar% {
    hotspotName = sprintf("%s_%s_%s", hotspotSigGenesDf[i, "space"], hotspotSigGenesDf[i, "start"], hotspotSigGenesDf[i, "end"])
    cat(sprintf("Processing %d of %d ...\n", i, nrow(hotspotSigGenesDf)))
    geneName = hotspotSigGenesDf[i, "SYMBOL"]
    cnvByGeneGeneDf = cnvByGeneDf[geneName,]
    if (!all(is.na(cnvByGeneGeneDf)) & nrow(cnvByGeneGeneDf) > 0) {
        cpids = colnames(cnvByGeneGeneDf[, cnvByGeneGeneDf > 0.3 | cnvByGeneGeneDf < -0.3])
    } else {
        cpids = c()
    }
    rsemLog2FoldGeneDf = subset(rsemLog2FoldPanDf, symbol == geneName)
    if (nrow(rsemLog2FoldGeneDf) > 0) {
        mpids = as.vector(sapply(strsplit(hotspotSigGenesDf[i, "sids"], ",", fixed = T)[[1]],
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
        hotspotSigGenesDf[i, "ttPvalueGeneLog2FoldExpMutVsWdt"] = ttPvalueGeneLog2FoldExpMutVsWdt
        hotspotSigGenesDf[i, "wtPvalueGeneLog2FoldExpMutVsWdt"] = wtPvalueGeneLog2FoldExpMutVsWdt
        if (!is.na(ttPvalueGeneLog2FoldExpMutVsWdt) & !is.na(wtPvalueGeneLog2FoldExpMutVsWdt)) {
            if (ttPvalueGeneLog2FoldExpMutVsWdt <= 0.05 | wtPvalueGeneLog2FoldExpMutVsWdt <= 0.05) {
                cat(sprintf("Found a significant expression change of %s in %s mutant group!\n", geneName, hotspotName))
                geneLog2FoldExpMutVsWdtDf = data.frame()
                geneLog2FoldExpMutVsWdtDf = rbind(geneLog2FoldExpMutVsWdtDf, data.frame(type=rep(hotspotName, length(mepids)), log2FoldRsem = rsemLog2FoldGeneMt))
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
        hotspotSigGenesDf[i, "ttPvalueGeneLog2FoldExpMutVsWdt"] = NA
        hotspotSigGenesDf[i, "wtPvalueGeneLog2FoldExpMutVsWdt"] = NA
    }
    hotspotSigGenesDf[i,]
}

head(hotspotSigGenesNewDf, 100)
#tail(hotspotSigGenesNewDf)

##
## Add annotation for known cancer genes
##
allOncoFile = file.path("/n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/allonco/allonco_20130923.tsv")
allOncoDf = read.delim(allOncoFile, header = T, as.is = T)
hotspotSigGenesNewDf$allOnco = ifelse(hotspotSigGenesNewDf$SYMBOL %in% allOncoDf$symbol, "Y", "N")

cgcFile = file.path("/n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/cosmic/cancer_gene_census.csv")
cgcDf = read.csv(cgcFile, header = T, as.is = T)
hotspotSigGenesNewDf$cosmic = ifelse(hotspotSigGenesNewDf$SYMBOL %in% cgcDf$Gene.Symbol, "Y", "N")

tcgaSmgFile = file.path("/n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/tcga-pancan/hcd.txt")
tcgaSmgDf = read.delim(tcgaSmgFile, header = T, as.is = T)
hotspotSigGenesNewDf$TCGA_PanCan_SMG = ifelse(hotspotSigGenesNewDf$SYMBOL %in% tcgaSmgDf$Gene.Symbol, "Y", "N")


##
## Add supporing evidence from HiC data?
##


##
## Save results to a file
##

hotspotSigVepOutMergedNewDf = merge(hotspotSigVepOutMergedDf, hotspotSigGenesNewDf, by = c("space", "start", "end", "SYMBOL", "sids"), all.x = TRUE)
hotspotSigVepOutMergedNewDf = hotspotSigVepOutMergedNewDf[order(-hotspotSigVepOutMergedNewDf$nrscnt),]
head(hotspotSigVepOutMergedNewDf)
hotspotSigVepOutMergedNewFile = file.path(panVcfDir, sprintf("TCGA_16_Cancer_Types.wgs.somatic.sanitized.hotspot%d.sig%d.vep_out.exp.with_scnt.new.txt", distCut, rscntPerSloci))
write.table(hotspotSigVepOutMergedNewDf, hotspotSigVepOutMergedNewFile, row.names = F, col.names = T, sep = "\t", quote = F)


##
save.image()
