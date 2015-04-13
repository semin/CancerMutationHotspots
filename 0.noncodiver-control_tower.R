rm(list=ls())

##
## Load libraries
##
require(doMC)
registerDoMC(10)

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
hotspotDir = file.path(vcfDir, "pancan")
hotspotDir = file.path(baseDir, "hotspot")
gdacDir = file.path(baseDir, "gdac")
gdacStdDataDir = file.path(gdacDir, "stddata__2015_02_04")
gdacAnlDataDir = file.path(gdacDir, "analyses__2014_10_17")
figDir = file.path(baseDir, "figure")
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

head(hotspotDf)
summary(hotspotDf$width)
hist(hotspotDf$width)
summary(hotspotDf$all_adj_scnt_p_value)
hist(hotspotDf$all_adj_scnt_p_value)

#hotspotFile = file.path(hotspotDir, sprintf("TCGA_16_Cancer_Types.wgs.somatic.sanitized.scnt.hotspot%d.txt", distCut))
#write.table(hotspotDf, hotspotFile, row.names = F, col.names = T, sep = "\t", quote = F)


##
## Filter significant hotspots
##

hotspotSigDf = subset(hotspotDf, all_adj_scnt_p_value < 0.05)
nrow(hotspotSigDf)
summary(hotspotSigDf$width)
hist(hotspotSigDf$width)

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
    print(hotspotChrVepSigFile)
    hotspotChrVepSigOutFile = gsub("vep_in", "vep_out", hotspotChrVepSigFile)
    hotspotChrVepSigLsfOutFile = gsub("txt", "lsfout", hotspotChrVepSigOutFile)
    system(sprintf("
bsub -g /nd/hotspot/vep \\
    -q i2b2_12h -W 12:0 \\
    -n 10 -R \"span[hosts=1]\" \\
    -o %s \\
    perl /home/sl279/vep/variant_effect_predictor.pl \\
        --fork 10 \\
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
        --custom /n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/ucsc/encode/wgEncodeAwgSegmentationCombinedH1hesc.bed.gz,wgEncodeAwgSegmentationChromhmmH1hesc,bed,overlap,0 \\
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
    write.table(hotspotSigVepOutMergedDf, hotspotSigVepMergedOutFile, row.names = F, col.names = T, sep = "\t", quote = F)
}


##
## Test mRNA expression changes 
##

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


## Test significance of expression changes of nearby genes
# ruby 3.noncodiver-hotspot_impact_on_expression.rb

hotspotAnnotChrFiles = mixedsort(Sys.glob(file.path(hotspotDir, sprintf("*.hotspot%d.fdr0.05.annotated.exptested.txt", distCut))))
hotspotAnnotDf = data.frame()
for (hotspotAnnotChrFile in hotspotAnnotChrFiles) {
    cat(sprintf("Reading %s ...\n", hotspotAnnotChrFile))
    hotspotAnnotChrDf = read.delim(hotspotAnnotChrFile, header = T, as.is = T)
    hotspotAnnotDf = rbind(hotspotAnnotDf, hotspotAnnotChrDf)
}
hotspotAnnotDf = hotspotAnnotDf[order(hotspotAnnotDf$all_adj_scnt_p_value),]
hotspotAnnotDf = hotspotAnnotDf[hotspotAnnotDf$all_adj_scnt_p_value < 0.01,]


##
## Add annotation for known cancer genes
##
allOncoFile = file.path("/n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/allonco/allonco_20130923.tsv")
allOncoDf = read.delim(allOncoFile, header = T, as.is = T)
hotspotAnnotDf$allOnco = ifelse(hotspotAnnotDf$SYMBOL %in% allOncoDf$symbol, "Y", "N")

intogenDriverSignalFile = file.path("/n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/intogen/Mutational_drivers_signals.tsv")
intogenDriverSignalDf = read.delim(intogenDriverSignalFile, header = T, as.is = T, comment.char = "#")
hotspotAnnotDf$intogenDriverSignals = ifelse(hotspotAnnotDf$SYMBOL %in% intogenDriverSignalDf$geneHGNCsymbol, "Y", "N")

intogenDriverFile = file.path("/n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/intogen/Drivers_type_role.tsv")
intogenDriverDf = read.delim(intogenDriverFile, header = T, as.is = T, comment.char = "#")
hotspotAnnotDf$intogenDrivers = ifelse(hotspotAnnotDf$SYMBOL %in% intogenDriverDf$geneHGNCsymbol, "Y", "N")

cgcFile = file.path("/n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/cosmic/cancer_gene_census.csv")
cgcDf = read.csv(cgcFile, header = T, as.is = T)
hotspotAnnotDf$cosmic = ifelse(hotspotAnnotDf$SYMBOL %in% cgcDf$Gene.Symbol, "Y", "N")

tcgaSmgFile = file.path("/n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver/tcga-pancan/hcd.txt")
tcgaSmgDf = read.delim(tcgaSmgFile, header = T, as.is = T)
hotspotAnnotDf$TCGA_PanCan_SMG = ifelse(hotspotAnnotDf$SYMBOL %in% tcgaSmgDf$Gene.Symbol, "Y", "N")


##
## Add supporing evidence from HiC-like data
##


##
## Save results!
##

hotspotAnnotDf$X4dGenomeHomoSapiens = NULL
hotspotAnnotFile = file.path(hotspotDir, sprintf("TCGA_16_Cancer_Types.wgs.somatic.sanitized.hotspot%d.fdr0.01.annotated.txt", distCut))
write.table(hotspotAnnotDf, hotspotAnnotFile, row.names = F, col.names = T, sep = "\t", quote = F)
save.image(imageFile)
