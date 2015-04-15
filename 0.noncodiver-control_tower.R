rm(list=ls())

##
## Load libraries
##
require(doMC)
registerDoMC(10)

require(ape)
require(Hmisc)
require(proxy)
require(sqldf)
require(gdata)
require(gtools)
require(gplots)
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
circosDir = file.path(baseDir, "circos")
imageFile = file.path(baseDir, "script/0.noncodiver-control_tower.RData")

chrs = c(1:22, "X", "Y")
chrs = factor(chrs, level=chrs)
cchrs = sapply(chrs, function(x) paste("chr", x, sep=""))
cchrs = factor(cchrs, level=cchrs)
distCut = 100

# Load chrom size information
hg19File = file.path(baseDir, "ucsc/database/hg19.genome")
hg19Df = read.delim(hg19File, header=T, as.is=T)
hg19Df = hg19Df[hg19Df$chrom %in% chrs,]
hg19Df$chrom = factor(hg19Df$chrom, levels=chrs)
hg19Df = hg19Df[order(hg19Df$chrom),]

##
## Detect hotspot candidates
##

# ruby 1.noncodiver-hotspot_detection.rb

##
## Merge hotspot candidates
##
hotspotChrFiles = mixedsort(Sys.glob(file.path(hotspotDir, sprintf("*.hotspot%d.txt", distCut))))
hotspotDf = data.frame()
for (hotspotChrFile in hotspotChrFiles) {
    cat(sprintf("Merging %s ...\n", hotspotChrFile))
    hotspotChrDf = read.delim(hotspotChrFile, header = T, as.is = T)
    hotspotDf = rbind(hotspotDf, hotspotChrDf)
}

hotspotDf$all_adj_scnt_p_value = p.adjust(hotspotDf$scnt_p_value, method = "fdr")
hotspotDf = hotspotDf[order(hotspotDf$all_adj_scnt_p_value),]

#head(hotspotDf)
#summary(hotspotDf$width)
#hist(hotspotDf$width)
#summary(hotspotDf$all_adj_scnt_p_value)
#hist(hotspotDf$all_adj_scnt_p_value)

#hotspotFile = file.path(hotspotDir, sprintf("TCGA_16_Cancer_Types.wgs.somatic.sanitized.scnt.hotspot%d.txt", distCut))
#write.table(hotspotDf, hotspotFile, row.names = F, col.names = T, sep = "\t", quote = F)


##
## Filter significant hotspots
##

hotspotSigDf = subset(hotspotDf, all_adj_scnt_p_value < 0.05)
#nrow(hotspotSigDf)
#summary(hotspotSigDf$all_adj_scnt_p_value)
#hist(hotspotSigDf$width)

for (chr in chrs) {
    chr = as.character(chr)
    print(chr)
    hotspotSigChrDf = subset(hotspotSigDf, space == chr)
    hotspotSigChrFile = file.path(hotspotDir, sprintf("TCGA_16_Cancer_Types.wgs.somatic.chr%s.sanitized.scnt.hotspot%d.fdr0.05.txt", chr, distCut))
    hotspotSigChrVepFile = file.path(hotspotDir, sprintf("TCGA_16_Cancer_Types.wgs.somatic.chr%s.sanitized.scnt.hotspot%d.fdr0.05.vep_in.txt", chr, distCut))
    write.table(hotspotSigChrDf, hotspotSigChrFile, row.names = F, col.names = T, sep = "\t", quote = F)
    write.table(hotspotSigChrDf[, c("space", "start", "end", "allele", "strand")],
                hotspotSigChrVepFile, row.names = F, col.names = F, sep = "\t", quote = F)
}


##
## Annotate significant hotspots with VEP
##
hotspotChrVepSigFiles = mixedsort(Sys.glob(file.path(hotspotDir, sprintf("*hotspot%d.fdr0.05.vep_in.txt", distCut))))
for (hotspotChrVepSigFile in hotspotChrVepSigFiles) {
    hotspotChrVepSigOutFile = gsub("vep_in", "vep_out", hotspotChrVepSigFile)
    hotspotChrVepSigLsfOutFile = gsub("txt", "lsfout", hotspotChrVepSigOutFile)
    if (file.exists(hotspotChrVepSigLsfOutFile)) next
    system(sprintf("
bsub -g /nd/hotspot/vep \\
    -q i2b2_12h -W 12:0 \\
    -n 4 -R \"span[hosts=1]\" \\
    -o %s \\
    perl /home/sl279/vep/variant_effect_predictor.pl \\
        --fork 4 \\
        --force_overwrite \\
        --offline \\
        --no_stats \\
        --everything \\
        --check_existing \\
        --total_length \\
        --allele_number \\
        --no_escape \\
        --per_gene \\
        --symbol \\
        --gencode_basic \\
        --assembly GRCh37 \\
        --dir /home/sl279/.vep \\
        --plugin UpDownDistance,1000000,1000000 \\
        --fasta /home/sl279/.vep/homo_sapiens/78_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/encode/distalDhsToPromoterDhs.bed.gz,distalDhsToPromoterDhs,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/encode/dhsToGeneExpression.bed.gz,dhsToGeneExpression,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/fantom5/fantom5EnhancerTssAssociations.bed.gz,fantom5EnhancerTssAssociations,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/dbsuper/all/allDbSuperEnhancerGeneAssociations.bed.gz,allDbSuperEnhancerGeneAssociations,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/4dgenome/4dGenomeHomoSapiens.bed.gz,fourDGenomeHomoSapiens,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/insitu-hic/insitu_HiC_GM12878_100kb_MAPQGE30_100kb_SQRTVC.bed.gz,insitu_HiC_GM12878_100kb_MAPQGE30_100kb_SQRTVC,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/insitu-hic/insitu_HiC_HMEC_100kb_intra_MAPQGE30_100kb_SQRTVC.bed.gz,insitu_HiC_HMEC_100kb_intra_MAPQGE30_100kb_SQRTVC,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/insitu-hic/insitu_HiC_HUVEC_100kb_intra_MAPQGE30_100kb_SQRTVC.bed.gz,insitu_HiC_HUVEC_100kb_intra_MAPQGE30_100kb_SQRTVC,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/insitu-hic/insitu_HiC_IMR90_100kb_intra_MAPQGE30_100kb_SQRTVC.bed.gz,insitu_HiC_IMR90_100kb_intra_MAPQGE30_100kb_SQRTVC,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/insitu-hic/insitu_HiC_K562_100kb_intra_MAPQGE30_100kb_SQRTVC.bed.gz,insitu_HiC_K562_100kb_intra_MAPQGE30_100kb_SQRTVC,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/insitu-hic/insitu_HiC_KBM7_100kb_intra_MAPQGE30_100kb_SQRTVC.bed.gz,insitu_HiC_KBM7_100kb_intra_MAPQGE30_100kb_SQRTVC,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/insitu-hic/insitu_HiC_NHEK_100kb_intra_MAPQGE30_100kb_SQRTVC.bed.gz,insitu_HiC_NHEK_100kb_intra_MAPQGE30_100kb_SQRTVC,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/encode/wgEncodeAwgDnaseMasterSites.bed.gz,wgEncodeAwgDnaseMasterSites,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/encode/wgEncodeRegDnaseClusteredV3.bed.gz,wgEncodeRegDnaseClusteredV3,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/fantom5/fantom5PermissiveEnhancers.bed.gz,fantom5PermissiveEnhancers,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/dbsuper/all/allDbSuperEnhancers.bed.gz,allDbSuperEnhancers,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/roadmap/dnase/roadmapDnaseDyadic.bed.gz,roadmapDnaseDyadic,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/roadmap/dnase/roadmapDnaseEnh.bed.gz,roadmapDnaseEnh,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/roadmap/dnase/roadmapDnaseProm.bed.gz,roadmapDnaseProm,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/encode/wgEncodeRegTfbsClusteredWithCellsV3.bed.gz,wgEncodeRegTfbsClusteredWithCellsV3,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/encode/wgEncodeAwgSegmentationCombinedGm12878.bed.gz,wgEncodeAwgSegmentationCombinedGm12878,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/encode/wgEncodeAwgSegmentationCombinedH1hesc.bed.gz,wgEncodeAwgSegmentationCombinedH1hesc,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/encode/wgEncodeAwgSegmentationCombinedHelas3.bed.gz,wgEncodeAwgSegmentationCombinedHelas3,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/encode/wgEncodeAwgSegmentationCombinedHepg2.bed.gz,wgEncodeAwgSegmentationCombinedHepg2,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/encode/wgEncodeAwgSegmentationCombinedHuvec.bed.gz,wgEncodeAwgSegmentationCombinedHuvec,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/encode/wgEncodeAwgSegmentationCombinedK562.bed.gz,wgEncodeAwgSegmentationCombinedK562,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/database/all_est.bed.gz,all_est,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/database/gwasCatalog.bed.gz,gwasCatalog,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/database/targetScanS.bed.gz,targetScanS,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/database/tfbsConsSites.bed.gz,tfbsConsSites,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/database/wgRna.bed.gz,wgRna,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/database/phastConsElements46way.bed.gz,phastConsElements46way,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/database/phastConsElements100way.bed.gz,phastConsElements100way,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/database/gap.bed.gz,gap,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/database/genomicSuperDups.bed.gz,genomicSuperDups,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/database/rmsk.bed.gz,rmsk,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/database/nestedRepeats.bed.gz,nestedRepeats,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/database/simpleRepeat.bed.gz,simpleRepeat,bed,overlap,0 \\
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/database/microsat.bed.gz,microsat,bed,overlap,0 \\
        --input_file %s \\
        --output_file %s
", hotspotChrVepSigLsfOutFile, hotspotChrVepSigFile, hotspotChrVepSigOutFile), intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE, wait = F)
}


##
## Expand vep_out files using Ruby script
##
# ➜  hotspot>  for f in *.vep_out.txt;do ruby ../script/expand_vep_out.rb $f > $f.expanded;done


##
## Annotate significant hotspots with VEP output
##
hotspotSigVepOutChrFiles = mixedsort(Sys.glob(file.path(hotspotDir, sprintf("TCGA_16_Cancer_Types.wgs.somatic.chr*.hotspot%d.fdr0.05.vep_out.txt.expanded", distCut, rscntPerSloci))))
for (hotspotSigVepOutChrFile in hotspotSigVepOutChrFiles) {
    cat(sprintf("Reading %s ...\n", hotspotSigVepOutChrFile))
    hotspotSigVepOutChrDf = read.delim(hotspotSigVepOutChrFile, header = T, as.is = T)
    hotspotSigVepOutChrDf = subset(hotspotSigVepOutChrDf, genomicSuperDups == "" & !grepl("trf", simpleRepeat))
    hotspotSigVepOutChrDf[, "space"] = gsub("(^\\S+?):.*", "\\1", hotspotSigVepOutChrDf$Location, perl = T)
    hotspotSigVepOutChrDf[, "start"] = sapply(hotspotSigVepOutChrDf$Location, function(x) { strsplit(strsplit(x, ":", fixed = T)[[1]][2], "-", fixed = T)[[1]][1] })
    hotspotSigVepOutChrDf[, "end"] = sapply(hotspotSigVepOutChrDf$Location, function(x) { strsplit(strsplit(x, ":", fixed = T)[[1]][2], "-", fixed = T)[[1]][2] })
    hotspotSigVepOutMergedDf = merge(hotspotSigVepOutChrDf, hotspotSigDf, by = c("space", "start", "end"), all.x = TRUE)
    hotspotSigVepOutMergedDf = hotspotSigVepOutMergedDf[order(hotspotSigVepOutMergedDf$all_adj_scnt_p_value),]
    hotspotSigVepMergedOutFile = gsub("vep_out.txt.expanded", "annotated.txt", hotspotSigVepOutChrFile, fixed = T)
    cat(sprintf("Writing %s ...\n", hotspotSigVepMergedOutFile))
    write.table(hotspotSigVepOutMergedDf, hotspotSigVepMergedOutFile, row.names = F, col.names = T, sep = "\t", quote = F)
}


##
## Test asssociation between nearby genes
##


## Load exome mutation data from GDAC
mafFiles = Sys.glob(file.path(gdacAnlDataDir, "*/*/*/*.maf"))
mafPanDf = data.frame()
for (i in 1:length(mafFiles)) {
    cat(sprintf("Reading %s ...\n", mafFile))
    mafDf = read.delim(mafFile, header = T, as.is = T)[, c(1:19)]
    if (i == 1) {
        mafColNames = colnames(mafDf)
    } else {
        colnames(mafDf) = mafColNames
    }
    mafDf$Cancer_Type = strsplit(basename(mafFile), "-")[[1]][1]
    mafPanDf = rbind(mafPanDf, mafDf)
}
nonFuncVarClasses = c("Silent")
mafPanFuncDf = subset(mafPanDf, Variant_Classification %nin% nonFuncVarClasses)
mafRdata = file.path(gdacAnlDataDir, "mafPanFunc.RData")
save(mafPanFuncDf, file = mafRdata)


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
rsemLog2FoldPanRdata = file.path(gdacStdDataDir, "rsemLog2FoldPan.RData")
save(rsemLog2FoldPanDf, file = rsemLog2FoldPanRdata)


## Load SCNA data for filter
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
cnvByGeneRdata = file.path(gdacStdDataDir, "cnvByGene.RData")
save(cnvByGeneDf, file = cnvByGeneRdata)


## Run ruby script for batch jobs
# ruby 2.noncodiver-hotspot_association.rb

## Collect association results
hotspotAnnotChrFiles = mixedsort(Sys.glob(file.path(hotspotDir, sprintf("*chr*.hotspot%d.fdr0.05.annotated.asstested.txt", distCut))))
#hotspotAnnotChrFiles = mixedsort(Sys.glob(file.path(hotspotDir, sprintf("*chr*.hotspot%d.fdr0.05.annotated.txt", distCut))))
hotspotAnnotDf = data.frame()
for (hotspotAnnotChrFile in hotspotAnnotChrFiles) {
    cat(sprintf("Reading %s ...\n", hotspotAnnotChrFile))
    hotspotAnnotChrDf = read.delim(hotspotAnnotChrFile, header = T, as.is = T)
    hotspotAnnotDf = rbind(hotspotAnnotDf, hotspotAnnotChrDf)
}
hotspotAnnotDf = hotspotAnnotDf[order(hotspotAnnotDf$all_adj_scnt_p_value),]
hotspotAnnotSigDf = hotspotAnnotDf[hotspotAnnotDf$all_adj_scnt_p_value < 0.01,]


##
## Add cancer distribution information
##
vcfFiles = Sys.glob(file.path(vcfDir, "cancer/*/*.stage4.vcf.gz"))
sampleIdsByCancerLst = list()
sampleIdToCancerDf = data.frame()
for (vcfFile in vcfFiles) {
    canerType = basename(dirname(vcfFile))
    sampleIds = samples(scanVcfHeader(vcfFile))
    sampleIdsByCancerLst[[canerType]] = sampleIds
    for (sampleId in sampleIds) {
        sampleIdToCancerDf = rbind(sampleIdToCancerDf, data.frame(id = sampleId, cancer = cancerType))
    }
}

hotspotAnnotSigDf$allele = NULL
hotspotAnnotSigDf$strand = NULL
hotspotAnnotSigGrpDf = sqldf('SELECT space, start, end, sids, all_adj_scnt_p_value FROM hotspotAnnotSigDf GROUP BY space, start, end ORDER BY space, start, end')

for (i in 1:nrow(hotspotAnnotSigGrpDf)) {
    sampleIds = strsplit(hotspotAnnotSigGrpDf[i, "sids"], ",")[[1]]
    for (c in names(sampleIdsByCancerLst)) {
        hotspotAnnotSigGrpDf[i, c] = sum(sampleIds %in% sampleIdsByCancerLst[[c]])
    }
}

hotspotAnnotSigGrpDf = hotspotAnnotSigGrpDf[order(hotspotAnnotSigGrpDf$all_adj_scnt_p_value),]
head(hotspotAnnotSigGrpDf)

sampleHotspotMutProfileDf <- foreach(sampleId=unlist(sampleIdsByCancerLst), .combine = rbind) %dopar% {
    print(sampleId)
    sapply(1:nrow(hotspotAnnotSigGrpDf), function(x) as.integer(grepl(sampleId, hotspotAnnotSigGrpDf[x,]$sids)))
}

head(sampleHotspotMutProfileDf)
rownames(sampleHotspotMutProfileDf) = unlist(sampleIdsByCancerLst)
sampleHotspotMutProfileDist = dist(sampleHotspotMutProfileDf)
sampleHotspotMutProfileHc = hclust(sampleHotspotMutProfileDist)
#plot(sampleHotspotMutProfileHc)

hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x, method="euclidean")

heatmap.2(as.matrix(sampleHotspotMutProfileDf),
          hclustfun=hclustfunc, distfun=distfunc,
          col = c("#FFFFFF", "#FF0000"),
          scale="none",
          trace="none",
          density.info="none",
          dendrogram = "both",
          labRow = "", labCol = "",
          Rowv = TRUE, Colv = TRUE)

#hotspotCancerProfileDf = hotspotAnnotSigGrpDf[, names(sampleIdsByCancerLst)]
#hotspotCancerProfileDf.2 = t(hotspotCancerProfileDf)
#hotspotCancerProfileDf.2.dist = dist(hotspotCancerProfileDf.2)
#hotspotCancerProfileDf.2.hc = hclust(hotspotCancerProfileDf.2.dist, method="complete")
#plot(hotspotCancerProfileDf.2.hc)
#plot(as.phylo(hotspotCancerProfileDf.2.hc), type='fan', label.offset=0.1, no.margin=TRUE)

#hotspotAnnotSigNewDf = merge(hotspotAnnotSigDf, hotspotAnnotSigGrpDf, by = c("space", "start", "end", "sids"), all.x = T)
#hotspotAnnotSigNewDf = hotspotAnnotSigNewDf[order(hotspotAnnotSigNewDf$all_adj_scnt_p_value),]



##
## Add annotation for known cancer genes
##
allOncoFile = file.path("/n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/allonco/allonco_20130923.tsv")
allOncoDf = read.delim(allOncoFile, header = T, as.is = T)
hotspotAnnotSigDf$allOnco = ifelse(hotspotAnnotSigDf$SYMBOL %in% allOncoDf$symbol, "Y", "N")

intogenDriverSignalFile = file.path("/n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/intogen/Mutational_drivers_signals.tsv")
intogenDriverSignalDf = read.delim(intogenDriverSignalFile, header = T, as.is = T, comment.char = "#")
hotspotAnnotSigDf$intogenDriverSignals = ifelse(hotspotAnnotSigDf$SYMBOL %in% intogenDriverSignalDf$geneHGNCsymbol, "Y", "N")

intogenDriverFile = file.path("/n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/intogen/Drivers_type_role.tsv")
intogenDriverDf = read.delim(intogenDriverFile, header = T, as.is = T, comment.char = "#")
hotspotAnnotSigDf$intogenDrivers = ifelse(hotspotAnnotSigDf$SYMBOL %in% intogenDriverDf$geneHGNCsymbol, "Y", "N")

cgcFile = file.path("/n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/cosmic/cancer_gene_census.csv")
cgcDf = read.csv(cgcFile, header = T, as.is = T)
hotspotAnnotSigDf$cosmic = ifelse(hotspotAnnotSigDf$SYMBOL %in% cgcDf$Gene.Symbol, "Y", "N")

tcgaSmgFile = file.path("/n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/tcga-pancan/hcd.txt")
tcgaSmgDf = read.delim(tcgaSmgFile, header = T, as.is = T)
hotspotAnnotSigDf$TCGA_PanCan_SMG = ifelse(hotspotAnnotSigDf$SYMBOL %in% tcgaSmgDf$Gene.Symbol, "Y", "N")


##
## Add supporing evidence from HiC-like data
##


##
## Save results!
##
hotspotAnnotFile = file.path(hotspotDir, sprintf("TCGA_16_Cancer_Types.wgs.somatic.sanitized.hotspot%d.fdr0.01.annotated.txt", distCut))
write.table(hotspotAnnotSigDf, hotspotAnnotFile, row.names = F, col.names = T, sep = "\t", quote = F)
save.image(imageFile)
load(imageFile)


##
## Plot linear genome-wide distribution of hotspots
##
hotspotAnnotChrFiles = mixedsort(Sys.glob(file.path(hotspotDir, sprintf("*chr**.hotspot%d.fdr0.05.annotated.txt", distCut))))
hotspotAnnotGrpDf = data.frame()
for (hotspotAnnotChrFile in hotspotAnnotChrFiles) {
    cat(sprintf("Reading %s ...\n", hotspotAnnotChrFile))
    hotspotAnnotChrDf = read.delim(hotspotAnnotChrFile, header = T, as.is = T)
    hotspotAnnotChrDf$allele = NULL
    hotspotAnnotChrDf$strand = NULL
    hotspotAnnotChrGrpDf = sqldf('SELECT space, start, end, chr_adj_scnt_p_value, wgEncodeAwgSegmentationCombinedGm12878 AS state FROM hotspotAnnotChrDf GROUP BY space, start, end ORDER BY space, start, end')
    hotspotAnnotChrGrpDf$log10Qvalue = -log10(hotspotAnnotChrGrpDf$chr_adj_scnt_p_value)
    hotspotAnnotChrGrpDf$label = ifelse(hotspotAnnotChrGrpDf$log10Qvalue > 3,
                                        sprintf("chr%s:%d-%d:%s", hotspotAnnotChrGrpDf$space, hotspotAnnotChrGrpDf$start, hotspotAnnotChrGrpDf$end, hotspotAnnotChrGrpDf$state), "")
    hotspotAnnotGrpDf = rbind(hotspotAnnotGrpDf, hotspotAnnotChrGrpDf)
    baseFontSize = 15
    breakSize = 10^7 * 2
    yBy = ceiling(max(-log10(hotspotAnnotChrDf$chr_adj_scnt_p_value)) / 5)
    p = ggplot(hotspotAnnotChrGrpDf, aes(x = start, y = -log10(chr_adj_scnt_p_value))) + 
        geom_hline(yintercept = seq(1, max(-log10(hotspotAnnotChrDf$chr_adj_scnt_p_value) + yBy), by=yBy), size = 0.5, linetype="dashed", alpha=.4) +
        geom_point(colour = "red", size = 2) + 
        #geom_line(colour = "grey") + 
        geom_text(aes(x = start, y = -log10(chr_adj_scnt_p_value), label = label, angle = 30, hjust = -0.03, vjust = -0.03, size = 9)) +
        theme(axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
              axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
              axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", hjust=1),
              axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(color = "black", fill="white"),
              legend.title = element_text(size = baseFontSize - 2, face="plain", family="sans", hjust = 0),
              legend.text  = element_text(size = baseFontSize - 2, face="plain", family="sans"),
              legend.direction = "horizontal",
              legend.position = "none") +
        scale_x_continuous(name=paste("\nGenomic position on chromosome ", hotspotAnnotChrGrpDf$space[1], " (Mb)", sep=""),
                           breaks=seq(1, max(hotspotAnnotChrGrpDf$end), by=breakSize),
                           limits=c(1, max(hotspotAnnotChrGrpDf$end) + breakSize),
                           labels=as.integer(seq(0, max(hotspotAnnotChrGrpDf$end), by=breakSize)/10^6)) +
        scale_y_continuous(name = "Significance of hotspot recurrence (-log10(P)) \n",
                           breaks=seq(1, max(-log10(hotspotAnnotChrGrpDf$chr_adj_scnt_p_value) + yBy), by=yBy),
                           limits=c(1, max(-log10(hotspotAnnotChrGrpDf$chr_adj_scnt_p_value) + 2*yBy)))
    hotspotAnnotChrGrpFile = gsub(".txt", ".pdf", hotspotAnnotChrFile, fixed = T)
    ggsave(filename = hotspotAnnotChrGrpFile, plot = p, width = 17, height = 6)
}


##
## Gernerate circos plot input data (hotspot p-values and labels)
##
circosLogPDf = hotspotAnnotGrpDf
circosLogPDf$logp = log10(circosLogPDf[,4])
circosLogPDf$space = paste("hs", circosLogPDf$space, sep = "")
circosLogPFile = file.path(circosDir, "hotspot.logp")
write.table(circosLogPDf[, c("space", "start", "end", "logp")], circosLogPFile, row.names = F, col.names = F, sep = "\t", quote = F)
circosLabelFile = file.path(circosDir, "hotspot.label")
write.table(subset(circosLogPDf, logp < -4)[, c("space", "start", "end", "label")], circosLabelFile, row.names = F, col.names = F, sep = "\t", quote = F)
