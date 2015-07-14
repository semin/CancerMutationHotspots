rm(list=ls())

##
## Load libraries
##
require(doMC)
registerDoMC(5)

require(ape)
require(vegan)
require(Hmisc)
require(proxy)
require(sqldf)
require(gdata)
require(scales)
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
## Custom functions
##
chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 

##
## Initialize global variables
##
rootDir = "/home/sl279"
baseDir = file.path(rootDir, "BiO/Research/Hotspot")
imageFile = file.path(baseDir, "script/0.noncodiver-control_tower.RData")
load(imageFile)
#save.image(imageFile)

vcfDir = file.path(baseDir, "vcf")
hotspotDir = file.path(baseDir, "hotspot")
sangDir = file.path(baseDir, "sanger")
icgcDir = file.path(baseDir, "icgc")
tcgaDir = file.path(baseDir, "tcga")
gdacDir = file.path(baseDir, "gdac")
gdacStdDataDir = file.path(gdacDir, "stddata__2015_02_04")
gdacAnlDataDir = file.path(gdacDir, "analyses__2014_10_17")
figDir = file.path(baseDir, "figure")
tableDir = file.path(baseDir, "table")
circosDir = file.path(baseDir, "circos")

chrs = c(1:22, "X")
chrs = factor(chrs, level=chrs)
cchrs = sapply(chrs, function(x) paste("chr", x, sep=""))
cchrs = factor(cchrs, level=cchrs)
distCut = 100


##
## Load chrom size information
##
hg19File = file.path(baseDir, "ucsc/database/hg19.genome")
hg19Df = read.delim(hg19File, header=T, as.is=T)
hg19Df = hg19Df[hg19Df$chrom %in% chrs,]
hg19Df$chrom = factor(hg19Df$chrom, levels=chrs)
hg19Df = hg19Df[order(hg19Df$chrom),]
hg19Df$start = 1
hg19Df$end = hg19Df$size
hg19Gr = with(hg19Df, GRanges(seqnames = Rle(chrom), IRanges(start = start, end = end), strand = "*"))


##
## Combine TCGA, Stratton, and ICGC mutation call sets
##

## Read ICGA release 18 WGS SNV call sets
icgcWgsSnvFiles = mixedsort(Sys.glob(file.path(baseDir, "icgc", "release_18", "*.wgs.snv.tsv")))
icgcWgsSnvAllDf = data.frame()
for (icgcWgsSnvFile in icgcWgsSnvFiles) {
    if (file.info(icgcWgsSnvFile)$size > 0) {
        cat(sprintf("Reading %s ...\n", icgcWgsSnvFile))
        icgcWgsSnvDf = read.delim(icgcWgsSnvFile, header = T, as.is = T)
        colnames(icgcWgsSnvDf) = c("icgc_mutation_id", "icgc_donor_id", "project_code", "icgc_specimen_id", "icgc_sample_id",
                                   "matched_icgc_sample_id", "submitted_sample_id", "submitted_matched_sample_id", 
                                   "chromosome", "chromosome_start", "chromosome_end", 
                                   "reference_genome_allele", "mutated_from_allele", "mutated_to_allele", 
                                   "consequence_type", "sequencing_strategy")
        icgcWgsSnvGrpDf = sqldf('SELECT * FROM icgcWgsSnvDf 
                                GROUP BY submitted_sample_id, chromosome, chromosome_start, chromosome_end,
                                reference_genome_allele, mutated_from_allele, mutated_to_allele')
        icgcWgsSnvAllDf = rbind(icgcWgsSnvAllDf, icgcWgsSnvGrpDf)
    }
}

icgcWgsSnvAllDf$cancer = gsub("(\\w+)-\\w+", "\\1", icgcWgsSnvAllDf$project_code)
icgcWgsSnvAllDf$source = "ICGC"

length(unique(icgcWgsSnvAllDf$submitted_sample_id))
length(unique(icgcWgsSnvAllDf$cancer))

## Read alt expanded TCGA mutation call sets
tcgaVcfFiles = Sys.glob(file.path(vcfDir, "cancer/*/*.stage4.vcf.gz"))
tcgaSampleIdsByCancerLst = list()
tcgaSampleIdToCancerDf = data.frame()
for (tcgaVcfFile in tcgaVcfFiles) {
    cat(sprintf("Reading %s ...\n", tcgaVcfFile))
    cancerType = basename(dirname(tcgaVcfFile))
    tcgaSampleIds = samples(scanVcfHeader(tcgaVcfFile))
    tcgaSampleIdsByCancerLst[[cancerType]] = tcgaSampleIds
    for (tcgaSampleId in tcgaSampleIds) {
        tcgaSampleIdToCancerDf = rbind(tcgaSampleIdToCancerDf, data.frame(sid = tcgaSampleId, cancer = cancerType))
    }
}
tcgaSampleIdToCancerDf$project = tcgaSampleIdToCancerDf$cancer

tcgaWgsSnvChrFiles = mixedsort(Sys.glob(file.path(tcgaDir, "*.snv.txt")))
tcgaWgsSnvAllDf = data.frame()
for (tcgaWgsSnvChrFile in tcgaWgsSnvChrFiles) {
    cat(sprintf("Reading %s ...\n", tcgaWgsSnvChrFile))
    tcgaWgsSnvChrDf = read.delim(tcgaWgsSnvChrFile, header = F, as.is = T)
    colnames(tcgaWgsSnvChrDf) = c("space", "pos", "ref", "alt", "sid")
    tcgaWgsSnvAllDf = rbind(tcgaWgsSnvAllDf, tcgaWgsSnvChrDf)
}
tcgaWgsSnvAllDf$sid = gsub("^H_\\w{2}", "TCGA", tcgaWgsSnvAllDf$sid)
tcgaWgsSnvAllDf$sid = gsub("345a06d6-fa5c-4674-a847-88a6b537cf3c", "TCGA-KN-8437-01A-11D-2310-10", tcgaWgsSnvAllDf$sid)
#grep("TCGA", tcgaWgsSnvAllDf$sid, invert = T, value = T)
tcgaWgsSnvAllDf$source = "TCGA"
tcgaWgsSnvAllDf = merge(tcgaWgsSnvAllDf, tcgaSampleIdToCancerDf, by = ("sid"), all.x = T)
#length(unique(tcgaWgsSnvAllDf$sid))

length(unique(tcgaWgsSnvAllDf$sid))
length(unique(tcgaWgsSnvAllDf$cancer))


## Read Stratton call sets
sangCanWgsSnvFiles = Sys.glob(file.path(sangDir, "somatic_mutation_data/*/*.wgs.txt"))
sangAllWgsSnvDf = data.frame()
for (sangCanWgsSnvFile in sangCanWgsSnvFiles) {
    if (file.info(sangCanWgsSnvFile)$size > 0) {
        cat(sprintf("Reading %s ...\n", sangCanWgsSnvFile))
        if (grepl("ALL", basename(sangCanWgsSnvFile))) { cancer = "ALL" }
        else if (grepl("AML", basename(sangCanWgsSnvFile))) { cancer = "ALL" }
        else if (grepl("Breast", basename(sangCanWgsSnvFile))) { cancer = "BRCA" }
        else if (grepl("CLL", basename(sangCanWgsSnvFile))) { cancer = "CLL" }
        else if (grepl("Liver", basename(sangCanWgsSnvFile))) { cancer = "LICA" }
        else if (grepl("Lung", basename(sangCanWgsSnvFile))) { cancer = "LUAD" }
        else if (grepl("B-cell", basename(sangCanWgsSnvFile))) { cancer = "BCL" }
        else if (grepl("Medullo", basename(sangCanWgsSnvFile))) { cancer = "MBL" }
        else if (grepl("Pancreas", basename(sangCanWgsSnvFile))) { cancer = "PACA" }
        else if (grepl("Pilocytic", basename(sangCanWgsSnvFile))) { cancer = "PA" }
        sangCanWgsSnvDf = read.delim(sangCanWgsSnvFile, header = F, as.is = T)[, c(3, 4, 6, 7, 1)]
        colnames(sangCanWgsSnvDf) = c("space", "pos", "ref", "alt", "sid")
        sangCanWgsSnvDf$cancer = cancer
        sangCanWgsSnvDf$project = cancer
        if (nrow(sangCanWgsSnvDf) > 0) {
            sangAllWgsSnvDf = rbind(sangAllWgsSnvDf, sangCanWgsSnvDf)
        }
    }
}
sangAllWgsSnvDf$source = "Sanger"
#length(unique(sangAllWgsSnvDf$sid))

length(unique(sangAllWgsSnvDf$sid))
length(unique(sangAllWgsSnvDf$cancer))


## cancer all SNVs into one data.frame
wgsSnvAllDf = data.frame()
wgsSnvAllDf = rbind(wgsSnvAllDf, data.frame(icgcWgsSnvAllDf[, c("source", "cancer", "project_code", "submitted_sample_id",
                                                                "chromosome", "chromosome_start", 
                                                                "reference_genome_allele", "mutated_to_allele")]))

colnames(wgsSnvAllDf) = c("source", "cancer", "project", "sid", "space", "pos", "ref", "alt")

wgsSnvAllDf = rbind(wgsSnvAllDf, tcgaWgsSnvAllDf[, colnames(wgsSnvAllDf)])
wgsSnvAllDf = rbind(wgsSnvAllDf, sangAllWgsSnvDf[, colnames(wgsSnvAllDf)])

nrow(wgsSnvAllDf)
length(unique(wgsSnvAllDf$sid))
length(unique(wgsSnvAllDf$cancer))


##
## Basic statistics for sample distribution and mutation rates
##
wgsSnvAllGrpByLocDf = sqldf('SELECT space, pos, ref, alt, COUNT(DISTINCT sid) AS scnt FROM wgsSnvAllDf GROUP BY space, pos, ref, alt')
nrow(wgsSnvAllGrpByLocDf)

wgsSnvAllGrpBySidDf = sqldf('SELECT source, cancer, project, sid, COUNT(*) mcnt FROM wgsSnvAllDf GROUP BY sid')
wgsSnvAllGrpBySidDf$mrate = 1000000 * wgsSnvAllGrpBySidDf$mcnt / 3036303846
wgsSnvAllGrpBySidDf = with(wgsSnvAllGrpBySidDf, wgsSnvAllGrpBySidDf[order(cancer, mcnt),])
wgsSnvAllGrpBySidGrpByCancerDf = sqldf('SELECT cancer, COUNT(sid) AS sid_cnt, AVG(mcnt) AS avg_mcnt, AVG(mrate) AS avg_mrate FROM wgsSnvAllGrpBySidDf GROUP BY cancer')

## Sample counts
baseFontSize = 11
wgsSnvAllGrpBySidGrpByCancerDf = with(wgsSnvAllGrpBySidGrpByCancerDf, wgsSnvAllGrpBySidGrpByCancerDf[order(cancer),])
p = ggplot(data = wgsSnvAllGrpBySidGrpByCancerDf, aes(x=factor(cancer, levels = wgsSnvAllGrpBySidGrpByCancerDf$cancer), y=sid_cnt)) +
    geom_bar(stat="identity") +
    theme(plot.title   = element_text(size = baseFontSize, face="bold"),
          axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
          axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
          axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 45, hjust=1),
          axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
          strip.text.x = element_text(size = baseFontSize, face="bold"),
          strip.text.y = element_text(size = baseFontSize, face="bold", angle = 0),
          legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
          legend.text  = element_text(size = baseFontSize, family="sans"),
          legend.direction = "horizontal",
          legend.position = "none") +
    scale_x_discrete(name="") +
    scale_y_continuous("# of samples\n")

wgsSnvAllGrpBySidGrpByCancerFile = file.path(baseDir, "figure", "Sample_count.wgs.pdf")
ggsave(filename = wgsSnvAllGrpBySidGrpByCancerFile, plot = p, width = 20, height = 5, unit = "cm")

## Mutation rates
wgsSnvAllGrpBySidGrpByCancerDf = with(wgsSnvAllGrpBySidGrpByCancerDf, wgsSnvAllGrpBySidGrpByCancerDf[order(avg_mcnt),])
baseFontSize = 11
p = ggplot(data = wgsSnvAllGrpBySidDf, aes(x=factor(cancer, levels = wgsSnvAllGrpBySidGrpByCancerDf$cancer), y=mrate)) +
    geom_boxplot() + 
    theme(plot.title   = element_text(size = baseFontSize, face="bold"),
          axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
          axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
          axis.text.x = element_text(size = baseFontSize, face="plain", family="sans", angle = 45, hjust = 1, colour = "black"),
          axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
          strip.text.x = element_text(size = baseFontSize, face="bold"),
          strip.text.y = element_text(size = baseFontSize, face="bold", angle = 0),
          legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
          legend.text  = element_text(size = baseFontSize - 5, family="sans"),
          legend.direction = "horizontal",
          legend.position = "bottom") +
    scale_x_discrete(name="") +
    scale_y_continuous("# of SNVs per Mb\n")

wgsSnvAllGrpBySidFile = file.path(baseDir, "figure", "Mutation_rates.wgs.pdf")
ggsave(filename = wgsSnvAllGrpBySidFile, plot = p, width = 30, height = 10, unit = "cm")



##
## Collapse combined SNVs
##
wgsSnvAllGrpDf = sqldf('SELECT space, pos, ref, alt, group_concat(sid) AS sids FROM wgsSnvAllDf GROUP BY space, pos, ref, alt')

for (chr in chrs) {
    wgsSnvAllGrpChrDf = subset(wgsSnvAllGrpDf, space == chr)
    wgsSnvAllGrpChrDf$scnt = sapply(wgsSnvAllGrpChrDf$sid, function(x) length(strsplit(x, ",")[[1]]))
    wgsSnvAllGrpChrFile = file.path(hotspotDir, sprintf("cancer.wgs.somatic.snv.chr%s.txt", chr))
    cat(sprintf("Writing %s ...\n", wgsSnvAllGrpChrFile))
    write.table(wgsSnvAllGrpChrDf, wgsSnvAllGrpChrFile, row.names = F, col.names = T, sep = "\t", quote = F)
}


##
## Detect hotspot candidates
##

# ruby 1.noncodiver-hotspot_detection.rb

##
## Merge hotspot candidates and calculate FDRs
##
hotspotChrFiles = mixedsort(Sys.glob(file.path(hotspotDir, sprintf("cancer.wgs.*.hotspot%d.txt", distCut))))
hotspotDf = data.frame()
for (hotspotChrFile in hotspotChrFiles) {
    cat(sprintf("Merging %s ...\n", hotspotChrFile))
    hotspotChrDf = read.delim(hotspotChrFile, header = T, as.is = T)
    hotspotDf = rbind(hotspotDf, hotspotChrDf)
}

subset(hotspotDf, is.na(scnt_p_value))

summary(hotspotDf$scnt_p_value)
summary(hotspotDf$mcnt_p_value)

hotspotDf$all_adj_scnt_p_value = p.adjust(hotspotDf$scnt_p_value, method = "fdr")
hotspotDf$all_adj_mcnt_p_value = p.adjust(hotspotDf$mcnt_p_value, method = "fdr")

hotspotDf = hotspotDf[order(hotspotDf$all_adj_mcnt_p_value),]
head(hotspotDf)
nrow(hotspotDf)


##
## Filter based on FDR (< 0.01) and generate chromosome-level hotspot files and VEP inputs files
##
hotspotSigDf = subset(hotspotDf, all_adj_mcnt_p_value < 0.01 & all_adj_scnt_p_value < 0.01)
nrow(hotspotSigDf)

for (chr in chrs) {
    chr = as.character(chr)
    print(chr)
    hotspotSigChrDf = subset(hotspotSigDf, space == chr)
    hotspotSigChrFile = file.path(hotspotDir, sprintf("cancer.wgs.somatic.chr%s.hotspot%d.fdr0.01.txt", chr, distCut))
    hotspotSigChrVepFile = file.path(hotspotDir, sprintf("cancer.wgs.somatic.chr%s.hotspot%d.fdr0.01.vep_in.txt", chr, distCut))
    write.table(hotspotSigChrDf, hotspotSigChrFile, row.names = F, col.names = T, sep = "\t", quote = F)
    write.table(hotspotSigChrDf[, c("space", "start", "end", "allele", "strand")], hotspotSigChrVepFile, row.names = F, col.names = F, sep = "\t", quote = F)
}


##
## Postprocess promoter capture results (http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2323/)
##
promIntFiles = Sys.glob(file.path(baseDir, "e-mtab-2323/*interactions.txt"))
for (promIntFile in promIntFiles) {
    cat(sprintf("Processing %s ...\n", promIntFile))
    promIntDf = read.delim(promIntFile, header = T, as.is = T)
    if (grepl("promoter-promoter", basename(promIntFile))) {
        promIntDf = read.delim(promIntFile, header = T, as.is = T)
        colnames(promIntDf) = NULL
        promIntDf1 = promIntDf[,c(1:3,10)]
        colnames(promIntDf1) = c("chr", "start", "end", "Symbol")
        promIntDf2 = promIntDf[,c(7:9,4)]
        colnames(promIntDf2) = c("chr", "start", "end", "Symbol")
        promIntBedDf = rbind(promIntDf1, promIntDf2)
    } else {
        promIntDf = read.delim(promIntFile, header = T, as.is = T)
        promIntBedDf = promIntDf[,c("chr", "start", "end", "Symbol")]
    }
    #for (i in 1:nrow(promIntBedDf)) {
        #promIntBedDf[i, "Symbol"] = paste(unique(sapply(strsplit(promIntBedDf[i, "Symbol"], "|", fixed = T)[[1]],
                                                        #function(x) gsub("(.*)-\\d{3}", "\\1", x))), collapse = "|")
    #}
    promIntBedDf$chr = factor(promIntBedDf$chr, levels = cchrs)
    promIntBedDf = with(promIntBedDf, promIntBedDf[order(chr, start, end),])
    promIntBedFile = gsub("txt", "bed", promIntFile)
    write.table(promIntBedDf, promIntBedFile, row.names = F, col.names = F, sep = "\t", quote = F)
}


##
## Postprocess insitu-HiC HiCCUPS loop lists
##
refGeneFile = file.path(baseDir, "ucsc/database/refGene.txt.gz")
refGeneDf = read.delim(gzfile(refGeneFile), header = F, as.is = T)
colnames(refGeneDf) = c("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd",
                        "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")
refGeneSelDf = refGeneDf[, c(2:6,13)]
colnames(refGeneSelDf) = c("transcript", "space", "strand", "start", "end", "gene")
refGeneSelDf$tss = with(refGeneSelDf, ifelse(strand == "+", start, end))
refGenePromDf = sqldf('SELECT gene, strand, space, tss FROM refGeneSelDf GROUP BY gene, tss')
prom_margin = 2000
refGenePromDf$prom_start = refGenePromDf$tss - prom_margin
refGenePromDf$prom_end = refGenePromDf$tss + prom_margin
refGenePromDf$space = gsub("chr", "", refGenePromDf$space)
refGenePromDf$space = factor(refGenePromDf$space, levels = chrs)
refGenePromDf = refGenePromDf[refGenePromDf$space %in% chrs,]
refGenePromGr = with(refGenePromDf, sort(GRanges(seqnames = Rle(space),
                                            IRanges(start = prom_start, end = prom_end),
                                            strand = strand,
                                            gene = gene,
                                            tss = tss)))

loopFiles = Sys.glob(file.path(baseDir, "insitu-hic", "*NHEK*looplist.txt"))
for (loopFile in loopFiles) {
    cat(sprintf("Processing %s ...\n", loopFile))
    loopDf = read.delim(loopFile, header = T, as.is = T)
    loopDf$dist = with(loopDf, y1 - x2)
    loopMirrDf = loopDf[, c("chr2", "y1", "y2", "chr1", "x1", "x2",
                            "color", "o", "e_bl", "e_donut", "e_h", "e_v",
                            "fdr_bl", "fdr_donut", "fdr_h", "fdr_v", "num_collapsed",
                            "centroid2", "centroid1", "radius", "dist")]
    colnames(loopMirrDf) = colnames(loopDf)
    loopCombDf = rbind(loopDf, loopMirrDf)
    loopCombDf$chr1 = factor(loopCombDf$chr1, levels = chrs)
    loopCombDf$chr2 = factor(loopCombDf$chr2, levels = chrs)
    targetDf = loopCombDf[, c("chr2", "y1", "y2")]
    targetNrDf = sqldf("SELECT * FROM targetDf GROUP BY chr2, y1, y2")
    targetNrGr = with(targetNrDf, sort(GRanges(seqnames = Rle(chr2),
                                               IRanges(start = y1, end = y2),
                                               strand = "*")))
    promHitsOls = findOverlaps(targetNrGr, refGenePromGr)
    for(i in unique(queryHits(promHitsOls))) {
        values(targetNrGr)[i, "genes"] = paste(unique(refGenePromGr[subjectHits(promHitsOls[queryHits(promHitsOls) == i])]$gene), collapse = "|")
    }
    targetNrDf = as.data.frame(targetNrGr)
    loopCombGenesDf = sqldf('SELECT l.*, t.genes 
                            FROM loopCombDf AS l 
                            LEFT JOIN targetNrDf AS t
                            ON l.chr2 == t.seqnames AND l.y1 == t.start AND l.y2 == t.end')
    loopCombGenesDf = with(loopCombGenesDf, loopCombGenesDf[order(chr1, x1, x2, chr2, y1, y2),])
    loopCombGenesFile = gsub(".txt", ".annotated.txt", loopFile)
    write.table(loopCombGenesDf, loopCombGenesFile, row.names = F, col.names = T, sep = "\t", quote = F)
    loopCombGenesBedFile = gsub(".txt", ".annotated.bed", loopFile)
    write.table(loopCombGenesDf[!is.na(loopCombGenesDf$genes), c("chr1", "x1", "x2", "genes")], loopCombGenesBedFile, row.names = F, col.names = F, sep = "\t", quote = F)
}


##
## Annotate hotspots with VEP
##
# ruby 2.noncodiver-hotspot_annotation.rb


##
## Expand vep_out files using Ruby script
##

# ➜  hotspot>  for f in *.vep_out.txt;do ruby ../script/2.noncodiver-expand_vep_out.rb $f > $f.expanded;done

##
## Read VEP output and filter out suspicious hotspots overlapping segmental duplications and simple repeats
##
hotspotSigVepOutChrFiles = mixedsort(Sys.glob(file.path(hotspotDir, sprintf("cancer*.chr*.hotspot%d.fdr0.01.vep_out.txt.expanded", distCut))))
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

rm(hotspotSigVepOutChrDf)
rm(hotspotSigVepOutMergedDf)


##
## Postprocess VEP annotated hotspots
##

hotspotAnnotChrFiles = mixedsort(Sys.glob(file.path(hotspotDir, sprintf("cancer*.chr*.hotspot%d.fdr0.01.annotated.txt", distCut))))
hotspotSigAnnotDf = data.frame()
hotspotSigAnnotGrpDf = data.frame()
for (hotspotAnnotChrFile in hotspotAnnotChrFiles) {
    cat(sprintf("Reading %s ...\n", hotspotAnnotChrFile))
    hotspotAnnotChrDf = read.delim(hotspotAnnotChrFile, header = T, as.is = T)
    hotspotSigAnnotDf = rbind(hotspotSigAnnotDf, hotspotAnnotChrDf)
    hotspotAnnotChrDf$allele = NULL
    hotspotAnnotChrDf$strand = NULL
    hotspotAnnotChrGrpDf = sqldf(sprintf("SELECT Location, all_adj_mcnt_p_value, all_adj_scnt_p_value, sids, %s
                                 FROM hotspotAnnotChrDf 
                                 GROUP BY Location ORDER BY all_adj_mcnt_p_value ASC", paste(grep("coreMarks", colnames(hotspotAnnotChrDf), value = T), collapse = ",")))
    hotspotSigAnnotGrpDf = rbind(hotspotSigAnnotGrpDf, hotspotAnnotChrGrpDf)
}

#hotspotSigAnnotDf = with(hotspotSigAnnotDf, hotspotSigAnnotDf[order(all_adj_mcnt_p_value),])
#nrow(hotspotSigAnnotGrpDf)


##
## Test asssociation between hotspots and nearby genes
##

## Load exome mutation data
mafFiles = Sys.glob(file.path(gdacAnlDataDir, "*/*/*/*.maf"))
mafPanDf = data.frame()
for (i in 1:length(mafFiles)) {
    mafFile = mafFiles[i]
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
mafPanFuncDf$sid = unlist(mclapply(mafPanFuncDf$Tumor_Sample_Barcode,
                                   function(x) {
                                       elems = strsplit(x, "-", fixed = T)[[1]][1:4]
                                       elems[4] = substring(elems[4], 1, 2)
                                       paste(elems, collapse = "_") }, mc.cores = numCores))

mafRdata = file.path(gdacAnlDataDir, "mafPanFunc.RData")
save(mafPanFuncDf, file = mafRdata)


## Load mRNA expression data
rsemFiles = Sys.glob(file.path(gdacStdDataDir, "*/*/*/*RSEM_all.txt"))
rpkmFiles = Sys.glob(file.path(gdacStdDataDir, "STAD/*/*/*RPKM.txt"))
geneExpFiles = c(rsemFiles, rpkmFiles)
geneExpPanDf = data.frame()
geneExpLog2FoldPanDf = data.frame()
for (i in 1:length(geneExpFiles)) {
    geneExpFile = geneExpFiles[i]
    cat(sprintf("Processing %s ...\n", geneExpFile))
    geneExpDf = read.delim(geneExpFile, header = T, as.is = T)
    colnames(geneExpDf)[1] = "id"
    colnames(geneExpDf)[2:ncol(geneExpDf)] = as.vector(sapply(colnames(geneExpDf)[2:ncol(geneExpDf)],
                                                              function(x) paste(strsplit(x, ".", fixed = T)[[1]][1:4], collapse = "_")))
    geneExpDf$symbol = sapply(geneExpDf$id, function(x) strsplit(x, "|", fixed = T)[[1]][1])
    geneExpDf$entrez = sapply(geneExpDf$id, function(x) gsub(".*?(\\d+).*", "\\1", strsplit(x, "|", fixed = T)[[1]][2]))
    geneExpLog2Df = log2(geneExpDf[, grep("TCGA", colnames(geneExpDf))] + 1)
    geneExpLog2FoldDf = geneExpLog2Df - rowMeans(geneExpLog2Df)
    geneExpLog2FoldDf = data.frame(id = geneExpDf$id, symbol = geneExpDf$symbol, entrez = geneExpDf$entrez, geneExpLog2FoldDf)
    subset(geneExpDf, symbol == "OR9I1")

    if (i == 1) {
        geneExpPanDf = geneExpDf
        geneExpLog2FoldPanDf = geneExpLog2FoldDf
    } else {
        geneExpDf$id = NULL
        geneExpLog2FoldDf$id = NULL
        geneExpPanDf = merge(geneExpPanDf, geneExpDf, by = c("symbol", "entrez"), all.x = TRUE)
        geneExpLog2FoldPanDf = merge(geneExpLog2FoldPanDf, geneExpLog2FoldDf, by = c("symbol", "entrez"), all.x = TRUE)
    }
}

geneExpLog2FoldPanRdata = file.path(gdacStdDataDir, "geneExpLog2FoldPan.RData")
save(geneExpLog2FoldPanDf, file = geneExpLog2FoldPanRdata)


## Load SCNA data 
cnvByGeneFiles = Sys.glob(file.path(gdacAnlDataDir, "*/*/*/all_data_by_genes.txt"))
cnvByGeneDf = data.frame()
for (i in 1:length(cnvByGeneFiles)) {
    cnvByGeneFile = cnvByGeneFiles[i]
    cat(sprintf("Reading %s ...\n", cnvByGeneFile))
    cnvByGeneTmpDf = read.delim(cnvByGeneFile, header = T, as.is = T)
    cnvByGeneTmpDf = cnvByGeneTmpDf[, c(-2,-3)]
    colnames(cnvByGeneTmpDf)[2:ncol(cnvByGeneTmpDf)] = as.vector(sapply(colnames(cnvByGeneTmpDf)[2:ncol(cnvByGeneTmpDf)],
                                                                        function(x) {
                                                                            elems = strsplit(x, ".", fixed = T)[[1]][1:4]
                                                                            elems[4] = substring(elems[4], 1, 2)
                                                                            paste(elems, collapse = "_")
                                                                        }))
    if (i == 1) {
        cnvByGeneDf = cnvByGeneTmpDf
    } else {
        cnvByGeneDf = cbind(cnvByGeneDf, cnvByGeneTmpDf[, grep("TCGA", colnames(cnvByGeneTmpDf))])
    }
}
rownames(cnvByGeneDf) = cnvByGeneDf$Gene.Symbol
cnvByGeneDf$Gene.Symbol = NULL
cnvByGeneRdata = file.path(gdacStdDataDir, "cnvByGene.RData")
save(cnvByGeneDf, file = cnvByGeneRdata)


## Split hotspot files into pieces for association tests
hotspotAnnotChrFiles = mixedsort(Sys.glob(file.path(hotspotDir, sprintf("*chr*.hotspot%d.fdr0.01.annotated.txt", distCut))))
for (hotspotAnnotChrFile in hotspotAnnotChrFiles) {
    cat(sprintf("Splitting %s ...\n", hotspotAnnotChrFile))
    hotspotAnnotChrDf = read.delim(hotspotAnnotChrFile, header = T, as.is = T)
    hotspotAnnotChrSplittedLst <- split(hotspotAnnotChrDf, hotspotAnnotChrDf$Location)
    for (i in 1:length(hotspotAnnotChrSplittedLst)) {
        hotspotAnnotChrIndDir = file.path(baseDir, "hotspot", "chrs", hotspotAnnotChrDf$space[1])
        dir.create(hotspotAnnotChrIndDir, recursive = TRUE, showWarnings = FALSE)
        hotspotAnnotChrIndFile = file.path(hotspotAnnotChrIndDir, gsub("txt", sprintf("%05d.txt", i), basename(hotspotAnnotChrFile)))
        cat(sprintf("Writing %s ...\n", hotspotAnnotChrIndFile))
        write.table(hotspotAnnotChrSplittedLst[[i]], hotspotAnnotChrIndFile, row.names = F, col.names = T, sep = "\t", quote = F)
    }
}


## Run ruby script for batch association jobs
# ruby 3.noncodiver-hotspot_association.rb


## Collect association results
hotspotAssFiles = mixedsort(Sys.glob(file.path(hotspotDir, sprintf("chrs/*/*chr*.hotspot%d.fdr0.01.annotated.*.ass.txt", distCut))))
hotspotAssDt <- foreach(hotspotAssFile=hotspotAssFiles, .combine=function(...) rbindlist(list(...)), .multicombine=TRUE) %dopar% {
    cat(sprintf("Reading %s ...\n", hotspotAssFile))
    if (file.info(hotspotAssFile)$size > 1) {
        hotspotAssDt = fread(hotspotAssFile, header = T)
    }
}

hotspotAssDf = as.data.frame(hotspotAssDt)
hotspotSigAssDf = subset(hotspotAssDf, log2FoldExpressionWilcoxPvalue < 0.01)
nrow(hotspotSigAssDf)
head(hotspotSigAssDf)


## Update hotspots with association results
roadmapChromStates = c("1_TssA", "2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", "8_ZNF/Rpts", "9_Het", "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk", "15_Quies", "Unk")
roadmapChromFillColors = c(rgb(255,0,0,maxColorValue=255), rgb(255,69,0,maxColorValue=255), rgb(50,205,50,maxColorValue=255), rgb(0,128,0,maxColorValue=255), rgb(0,100,0,maxColorValue=255), rgb(194,225,5,maxColorValue=255), rgb(255,255,0,maxColorValue=255), rgb(102,205,170,maxColorValue=255), rgb(138,145,208,maxColorValue=255), rgb(205,92,92,maxColorValue=255), rgb(233,150,122,maxColorValue=255), rgb(189,183,107,maxColorValue=255), rgb(128,128,128,maxColorValue=255), rgb(192,192,192,maxColorValue=255), rgb(255,255,255,maxColorValue=255), rgb(000,000,000,maxColorValue=255))
names(roadmapChromFillColors) = roadmapChromStates
roadmapChromBorderColors = c(rgb(255,0,0,maxColorValue=255), rgb(255,69,0,maxColorValue=255), rgb(50,205,50,maxColorValue=255), rgb(0,128,0,maxColorValue=255), rgb(0,100,0,maxColorValue=255), rgb(194,225,5,maxColorValue=255), rgb(255,255,0,maxColorValue=255), rgb(102,205,170,maxColorValue=255), rgb(138,145,208,maxColorValue=255), rgb(205,92,92,maxColorValue=255), rgb(233,150,122,maxColorValue=255), rgb(189,183,107,maxColorValue=255), rgb(128,128,128,maxColorValue=255), rgb(192,192,192,maxColorValue=255), rgb(0,0,0,maxColorValue=255), rgb(0,0,0,maxColorValue=255))
names(roadmapChromBorderColors) = roadmapChromStates
roadmapChromStateLabels = c("Active TSS", "Flanking Active TSS", "Transcr. at gene 5' and 3'", "Strong transcription", "Weak transcription", "Genic enhancers", "Enhancers", "ZNF genes & repeats", "Heterochromatin", "Bivalent/Poised TSS", "Flanking Bivalent TSS/Enh", "Bivalent Enhancer", "Repressed PolyComb", "Weak Repressed PolyComb", "Quiescent/Low", "Unknown")
names(roadmapChromStateLabels) = roadmapChromStates

locs = unique(hotspotSigAnnotDf$Location)
locChunks = chunk(1:length(locs), 10)
hotspotSigAssAnnotDf <- foreach(j=1:length(locChunks), .combine = rbind) %dopar% {
    cat(sprintf("Processing chunk %d out of %d ...\n", j, length(locChunks)))
    hotspotSigAssAnnotChunkDf = data.frame()
    for (i in locChunks[[j]]) {
        loc = locs[i]
        hotspotIndAnnotDf = subset(hotspotSigAnnotDf, Location == loc)
        olGenes = subset(hotspotIndAnnotDf, Feature_type == "Transcript")$SYMBOL
        cisGenes = subset(hotspotSigAssDf, Location == loc & supportedBy != "GSE63525_GM12878_100kb_inter_MAPQGE30_SQRTVC")$Gene
        transGenes = subset(hotspotSigAssDf, Location == loc & supportedBy == "GSE63525_GM12878_100kb_inter_MAPQGE30_SQRTVC")$Gene
        allStates = mixedsort(unique(unlist(sapply(hotspotIndAnnotDf[1, grep("coreMarks", colnames(hotspotIndAnnotDf))], function(x) strsplit(x, ",")[[1]]))))
        hotspotColor = roadmapChromFillColors[allStates[1]]
        hotspotState = paste(allStates, collapse = ",")

        hotspotSigAnnotSubDf = subset(hotspotSigAnnotDf, Location == loc)
        splice_site_cnt = sum(grepl("splice", hotspotSigAnnotSubDf$Consequence))
        exon_cnt = sum(grepl("frame", hotspotSigAnnotSubDf$Consequence) | hotspotSigAnnotSubDf$EXON != "")
        five_prime_utr_cnt = sum(grepl("5_prime", hotspotSigAnnotSubDf$Consequence))
        three_prime_utr_cnt = sum(grepl("3_prime", hotspotSigAnnotSubDf$Consequence))
        intron_cnt = sum(grepl("intron", hotspotSigAnnotSubDf$Consequence))
        tot_region_cnt = splice_site_cnt + exon_cnt + five_prime_utr_cnt + three_prime_utr_cnt + intron_cnt

        genomic_types = c()
        if (exon_cnt > 0) genomic_types = c(genomic_types, "exon")
        if (five_prime_utr_cnt > 0) genomic_types = c(genomic_types, "five_prime_utr")
        if (three_prime_utr_cnt > 0) genomic_types = c(genomic_types, "three_prime_utr")
        if (splice_site_cnt > 0) genomic_types = c(genomic_types, "splice_site") 
        if (intron_cnt > 0) genomic_types = c(genomic_types, "intron")
        if (tot_region_cnt == 0) genomic_types = c(genomic_types, "igr")
        genomic_type = paste(genomic_types, collapse = ",")

        hotspotSigAssAnnotChunkDf = with(hotspotIndAnnotDf[1,], 
                                         rbind(hotspotSigAssAnnotChunkDf,
                                               data.frame(Location = loc,
                                                          space = space,
                                                          start = start,
                                                          end = end,
                                                          mcnt = mcnt,
                                                          scnt = scnt,
                                                          sids = sids,
                                                          all_adj_mcnt_p_value = all_adj_mcnt_p_value,
                                                          all_adj_scnt_p_value = all_adj_scnt_p_value,
                                                          types = genomic_type,
                                                          overlapping_genes = paste(olGenes, collapse = ","),
                                                          cis_genes = paste(cisGenes, collapse = ","),
                                                          trans_genes = paste(transGenes, collapse = ","),
                                                          states = hotspotState,
                                                          rep_state = allStates[1],
                                                          color = hotspotColor,
                                                          stringsAsFactors = FALSE)))
    }
    hotspotSigAssAnnotChunkDf
}

delta = .Machine$double.xmin
hotspotSigAssAnnotDf$all_adj_mcnt_p_value_d = hotspotSigAssAnnotDf$all_adj_mcnt_p_value + delta
hotspotSigAssAnnotDf$all_adj_scnt_p_value_d = hotspotSigAssAnnotDf$all_adj_scnt_p_value + delta

sum(as.numeric(hotspotSigAssAnnotDf$mcnt))
hotspotSigAssAnnotDf$width = hotspotSigAssAnnotDf$end - hotspotSigAssAnnotDf$start + 1
sum(as.numeric(hotspotSigAssAnnotDf$width)) / sum(as.numeric(hg19Df$size))

##
## Plot linear genome-wide distribution of hotspots
##
for (chr in chrs) {
    cat(sprintf("Generating a plot for the linear distribution of hotspots in chr%s ...\n", chr))
    hotspotSigAssAnnotChrDf = subset(hotspotSigAssAnnotDf, space == chr)
    hotspotSigAssAnnotChrDf$log10Qvalue = -log10(hotspotSigAssAnnotChrDf$all_adj_mcnt_p_value_d)
    for (i in 1:nrow(hotspotSigAssAnnotChrDf)) {
        if (hotspotSigAssAnnotChrDf[i, "log10Qvalue"] > 7) {
            olGeneNames = strsplit(hotspotSigAssAnnotChrDf[i, "overlapping_genes"], ",")[[1]]
            cisGeneNames = strsplit(hotspotSigAssAnnotChrDf[i, "cis_genes"], ",")[[1]]
            transGeneNames = strsplit(hotspotSigAssAnnotChrDf[i, "trans_genes"], ",")[[1]]
            allGeneNames = unique(c(olGeneNames, cisGeneNames))
            geneCnt = length(allGeneNames)
            if (geneCnt > 0) {
                hotspotSigAssAnnotChrDf[i, "label"] = ifelse(geneCnt > 7, 
                                                        paste(paste(allGeneNames[1:ceiling(geneCnt/2)], collapse = ","),
                                                                paste(allGeneNames[(ceiling(geneCnt/2)+1):geneCnt], collapse = ","), sep = "\n"),
                                                    paste(allGeneNames, collapse = ",")) 
            } else {
                hotspotSigAssAnnotChrDf[i, "label"] = hotspotSigAssAnnotChrDf[i, "Location"]
            }
        } else {
            hotspotSigAssAnnotChrDf[i, "label"] = NA
        }
    }
    maxGeneCnt = max(sapply(hotspotSigAssAnnotChrDf$label, function(x) length(strsplit(x, ",")[[1]])))
    baseFontSize = 15
    breakSize = 10^7 * 2
    yBy = ceiling(max(-log10(hotspotSigAssAnnotChrDf$all_adj_mcnt_p_value_d)) / 5)
    yByScale = maxGeneCnt 
    p = ggplot(hotspotSigAssAnnotChrDf, aes(x = start, y = -log10(all_adj_mcnt_p_value_d))) + 
        geom_hline(yintercept = seq(1, max(-log10(hotspotSigAssAnnotChrDf$all_adj_mcnt_p_value_d) + yBy), by=yBy), size = 0.5, linetype="dashed", alpha=.4) +
        geom_point(colour = "red", size = 2) + 
        geom_text(aes(x = start, y = -log10(all_adj_mcnt_p_value_d), label = hotspotSigAssAnnotChrDf$label, 
                      colour = factor(hotspotSigAssAnnotChrDf$rep_state), angle = 45, hjust = -0.01, vjust = -0.01, size = 9)) +
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
        scale_x_continuous(name=paste("\nGenomic position on chromosome ", hotspotSigAssAnnotChrDf$space[1], " (Mb)", sep=""),
                           breaks=seq(1, max(hotspotSigAssAnnotChrDf$end), by=breakSize),
                           limits=c(1, max(hotspotSigAssAnnotChrDf$end) + breakSize),
                           labels=as.integer(seq(0, max(hotspotSigAssAnnotChrDf$end), by=breakSize)/10^6)) +
        scale_y_continuous(name = "Significance of hotspot mutation (-log10(q)) \n",
                           breaks=seq(1, max(-log10(hotspotSigAssAnnotChrDf$all_adj_mcnt_p_value_d) + yBy), by=yBy),
                           limits=c(1, max(-log10(hotspotSigAssAnnotChrDf$all_adj_mcnt_p_value_d) + yByScale * yBy))) +
        scale_colour_manual(values = roadmapChromFillColors)
    hotspotAnnotChrGrpFile = file.path(baseDir, "figure", sprintf("Linear_distribution_of_hotspot100_fdr0.01.chr%s.pdf", chr))
    ggsave(filename = hotspotAnnotChrGrpFile, plot = p, width = 17, height = 6)
}


##
## Gernerate circos plot input data (hotspot p-values and labels)
##
col2rgbLabel = function(cl) apply(col2rgb(sapply(cl, function(x)x[1])), 2, function(n) paste(n, collapse = ","))
circosLogPDf = hotspotSigAssAnnotDf
circosLogPDf$log10Qvalue = log10(circosLogPDf$all_adj_mcnt_p_value_d)
circosLogPDf = circosLogPDf[order(circosLogPDf$log10Qvalue),]
circosLogPDf$space = paste("hs", circosLogPDf$space, sep = "")
for (i in 1:nrow(circosLogPDf)) {
    if (circosLogPDf[i, "log10Qvalue"] < -5) {
        olGeneNames = strsplit(circosLogPDf[i, "overlapping_genes"], ",")[[1]]
        cisGeneNames = strsplit(circosLogPDf[i, "cis_genes"], ",")[[1]]
        allGeneNames = unique(c(olGeneNames, cisGeneNames))
        geneCnt = length(allGeneNames)
        if (geneCnt > 0) {
            circosLogPDf[i, "label"] = paste(allGeneNames[1], collapse = ",")
        } else {
            circosLogPDf[i, "label"] = circosLogPDf[i, "Location"]
        }
    } else {
        circosLogPDf[i, "label"] = ""
    }
    circosLogPDf[i, "rgb"] = sprintf("color=%s", col2rgbLabel(circosLogPDf[i, "color"]))
}
circosLogPFile = file.path(circosDir, "hotspot.logq")
write.table(circosLogPDf[, c("space", "start", "end", "log10Qvalue")], circosLogPFile, row.names = F, col.names = F, sep = "\t", quote = F)
circosLabelFile = file.path(circosDir, "hotspot.label")
#write.table(subset(circosLogPDf, log10Qvalue < -7)[, c("space", "start", "end", "label", "rgb")], circosLabelFile, row.names = F, col.names = F, sep = "\t", quote = F)
write.table(circosLogPDf[1:300, c("space", "start", "end", "label", "rgb")], circosLabelFile, row.names = F, col.names = F, sep = "\t", quote = F)


##
## Distribution and correlation between significance of mutation rates and recurrence level
##
colnames(hotspotSigAssAnnotDf)
hist(hotspotSigAssAnnotDf$all_adj_mcnt_p_value_d)
hist(hotspotSigAssAnnotDf$all_adj_scnt_p_value_d)
with(hotspotSigAssAnnotDf, plot(-log2(all_adj_mcnt_p_value_d), -log2(all_adj_scnt_p_value_d), xlim = c(0, 100), ylim = c(0, 100)))


##
## Check enrichment of hotspots in genomic and epigenomic features
##
source("http://www.pmc.ucsc.edu/~mclapham/Rtips/G%20test.txt")

# Genomic features
refGeneDf$space = gsub("chr", "", refGeneDf$chrom)
refGeneDf = refGeneDf[refGeneDf$space %in% chrs,]
refGeneDf$space = factor(refGeneDf$space, levels = chrs)
refGeneDf$FiveUtrStart = ifelse(refGeneDf$strand == "+", refGeneDf$txStart, refGeneDf$cdsEnd)
refGeneDf$FiveUtrEnd = ifelse(refGeneDf$strand == "+", refGeneDf$cdsStart, refGeneDf$txEnd)
refGeneDf$ThreeUtrStart = ifelse(refGeneDf$strand == "+", refGeneDf$cdsEnd, refGeneDf$txStart)
refGeneDf$ThreeUtrEnd = ifelse(refGeneDf$strand == "+", refGeneDf$txEnd, refGeneDf$cdsStart)

chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
refGeneChunks = chunk(1:nrow(refGeneDf), 10)
refGeneExonSpliceDf <- foreach(j=1:length(refGeneChunks), .combine = rbind) %dopar% {
    cat(sprintf("Processing RefGene chunk %d/%d ...\n", j, length(refGeneChunks)))
    refGeneExonSpliceChunkDf = data.frame()
    for (i in refGeneChunks[[j]]) {
        exonStartPos = as.numeric(strsplit(refGeneDf[i, "exonStarts"], ",", fixed = T)[[1]])
        exonEndPos = as.numeric(strsplit(refGeneDf[i, "exonEnds"], ",", fixed = T)[[1]])
        refGeneExonSpliceChunkDf = with(refGeneDf[i,], rbind(refGeneExonSpliceChunkDf,
                                                             data.frame(type = "exon", space = space, strand = strand, start = exonStartPos, end = exonEndPos)))
        if (refGeneDf[i, "exonCount"] > 1) {
            if (refGeneDf[i, "strand"] == "+") {
                exonDStarts = exonEndPos[1:(length(exonEndPos)-1)]
                exonAStarts = exonStartPos[2:length(exonStartPos)]
                spliceDStarts = exonDStarts - 2
                spliceDEnds = exonDStarts + 8
                spliceAStarts = exonAStarts - 8
                spliceAEnds = exonAStarts + 2
            } else if (refGeneDf[i, "strand"] == "-") {
                exonAStarts = exonEndPos[1:(length(exonEndPos)-1)]
                exonDStarts = exonStartPos[2:length(exonStartPos)]
                spliceDStarts = exonDStarts - 8
                spliceDEnds = exonDStarts + 2
                spliceAStarts = exonAStarts - 2
                spliceAEnds = exonAStarts + 8
            }
            refGeneExonSpliceChunkDf = with(refGeneDf[i,], rbind(refGeneExonSpliceChunkDf, 
                                                                 data.frame(type = "splice", space = space, strand = strand, start = spliceDStarts, end = spliceDEnds)))
            refGeneExonSpliceChunkDf = with(refGeneDf[i,], rbind(refGeneExonSpliceChunkDf, 
                                                                 data.frame(type = "splice", space = space, strand = strand, start = spliceAStarts, end = spliceAEnds)))
        }
    }
    refGeneExonSpliceChunkDf
}

refGene5utrGr = unique(sort(with(refGeneDf, GRanges(seqnames = space, IRanges(start = FiveUtrStart, end = FiveUtrEnd), strand = "*"))))
refGene3utrGr = unique(sort(with(refGeneDf, GRanges(seqnames = space, IRanges(start = ThreeUtrStart, end = ThreeUtrEnd), strand = "*"))))
refGeneExonDf = subset(refGeneExonSpliceDf, type == "exon")
refGeneExonGr = unique(sort(with(refGeneExonDf, GRanges(seqnames = Rle(space), IRanges(start = start, end = end), strand = "*"))))
refGeneSpliceDf = subset(refGeneExonSpliceDf, type == "splice")
refGeneSpliceGr = unique(sort(with(refGeneSpliceDf, GRanges(seqnames = Rle(space), IRanges(start = start, end = end), strand = "*"))))
refGeneGenicGr = unique(sort(with(refGeneDf, GRanges(seqnames = space, IRanges(start = txStart, end = txEnd), strand = "*"))))
refGeneIgrGr = setdiff(hg19Gr, refGeneGenicGr)
refGeneIntronGr = setdiff(refGeneGenicGr, refGeneExonGr)

hotspotAbbDf = hotspotSigAssAnnotDf[, c("space", "start", "end", "mcnt", "scnt", "all_adj_mcnt_p_value", "all_adj_scnt_p_value")]
hotspotAbbDf = hotspotAbbDf[hotspotAbbDf$space %in% chrs,]
hotspotAbbGr = sort(as(as(hotspotAbbDf, "RangedData"), "GRanges"))

hotspotAbb5UtrIntGr = intersect(hotspotAbbGr, refGene5utrGr)
hotspotAbb5UtrContMat = matrix(c(sum(as.numeric(width(hotspotAbb5UtrIntGr))),
                                sum(as.numeric(width(refGene5utrGr))) - sum(as.numeric(width(hotspotAbb5UtrIntGr))),
                                sum(as.numeric(width(hotspotAbbGr))) - sum(as.numeric(width(hotspotAbb5UtrIntGr))),
                                sum(as.numeric(width(hg19Gr))) - sum(as.numeric(width(hotspotAbbGr))) - sum(as.numeric(width(refGene5utrGr))) + sum(as.numeric(width(hotspotAbb5UtrIntGr)))),
                              nrow = 2, dimnames = list(Hotspot = c("Hotpot", "Non-hotspot"), Feature = c("5'UTR", "Non-5'UTR")))
hotspotAbb5UtrIntG = g.test(hotspotAbb5UtrContMat)

hotspotAbb5UtrIntDf = data.frame(hotspot = c("Hotspot", "Non-hotspot"),
                                 feature = c("5' UTR", "5' UTR"),
                                 size = c(sum(as.numeric(width(hotspotAbb5UtrIntGr))),
                                           sum(as.numeric(width(refGene5utrGr))) - sum(as.numeric(width(hotspotAbb5UtrIntGr)))),
                                 fraction = c(sum(as.numeric(width(hotspotAbb5UtrIntGr))) / sum(as.numeric(width(hotspotAbbGr))),
                                              (sum(as.numeric(width(refGene5utrGr))) - sum(as.numeric(width(hotspotAbb5UtrIntGr)))) / sum(as.numeric(width(hg19Gr)))))

hotspotAbb5UtrOrP = poisson.test(hotspotAbb5UtrContMat[,1], hotspotAbb5UtrContMat[,2], alternative = "two.sided")
hotspotAbb5UtrOrR = hotspotAbb5UtrContMat[,1] / hotspotAbb5UtrContMat[,2]
hotspotAbb5UtrOrDf = data.frame(feature = "5' UTR",
                                odds_ratio = hotspotAbb5UtrOrR[1] / hotspotAbb5UtrOrR[2],
                                poisson_p = hotspotAbb5UtrOrP$p.value)


hotspotAbb3UtrIntGr = intersect(hotspotAbbGr, refGene3utrGr)
hotspotAbb3UtrContMat = matrix(c(sum(as.numeric(width(hotspotAbb3UtrIntGr))),
                                sum(as.numeric(width(refGene3utrGr))) - sum(as.numeric(width(hotspotAbb3UtrIntGr))),
                                sum(as.numeric(width(hotspotAbbGr))) - sum(as.numeric(width(hotspotAbb3UtrIntGr))),
                                sum(as.numeric(width(hg19Gr))) - sum(as.numeric(width(hotspotAbbGr))) - sum(as.numeric(width(refGene3utrGr))) + sum(as.numeric(width(hotspotAbb3UtrIntGr)))),
                              nrow = 2, dimnames = list(Hotspot = c("Hotpot", "Non-hotspot"), Feature = c("3'UTR", "Non-3'UTR")))
hotspotAbb3UtrIntG = g.test(hotspotAbb3UtrContMat)

hotspotAbb3UtrIntDf = data.frame(hotspot = c("Hotspot", "Non-hotspot"),
                                 feature = c("3' UTR", "3' UTR"),
                                 size = c(sum(as.numeric(width(hotspotAbb3UtrIntGr))),
                                           sum(as.numeric(width(refGene3utrGr))) - sum(as.numeric(width(hotspotAbb3UtrIntGr)))),
                                 fraction = c(sum(as.numeric(width(hotspotAbb3UtrIntGr))) / sum(as.numeric(width(hotspotAbbGr))),
                                              (sum(as.numeric(width(refGene3utrGr))) - sum(as.numeric(width(hotspotAbb3UtrIntGr)))) / sum(as.numeric(width(hg19Gr)))))

hotspotAbb3UtrOrP = poisson.test(hotspotAbb3UtrContMat[,1], hotspotAbb3UtrContMat[,2], alternative = "two.sided")
hotspotAbb3UtrOrR = hotspotAbb3UtrContMat[,1] / hotspotAbb3UtrContMat[,2]
hotspotAbb3UtrOrDf = data.frame(feature = "3' UTR",
                                odds_ratio = hotspotAbb3UtrOrR[1] / hotspotAbb3UtrOrR[2],
                                poisson_p = hotspotAbb3UtrOrP$p.value)


hotspotAbbExonIntGr = intersect(hotspotAbbGr, refGeneExonGr)
hotspotAbbExonContMat = matrix(c(sum(as.numeric(width(hotspotAbbExonIntGr))),
                                sum(as.numeric(width(refGeneExonGr))) - sum(as.numeric(width(hotspotAbbExonIntGr))),
                                sum(as.numeric(width(hotspotAbbGr))) - sum(as.numeric(width(hotspotAbbExonIntGr))),
                                sum(as.numeric(width(hg19Gr))) - sum(as.numeric(width(hotspotAbbGr))) - sum(as.numeric(width(refGeneExonGr))) + sum(as.numeric(width(hotspotAbbExonIntGr)))),
                              nrow = 2, dimnames = list(Hotspot = c("Hotpot", "Non-hotspot"), Feature = c("Exon", "Non-exon")))

hotspotAbbExonIntG = g.test(hotspotAbbExonContMat)

hotspotAbbExonIntDf = data.frame(hotspot = c("Hotspot", "Non-hotspot"),
                                 feature = c("Exon", "Exon"),
                                 size = c(sum(as.numeric(width(hotspotAbbExonIntGr))),
                                           sum(as.numeric(width(refGeneExonGr))) - sum(as.numeric(width(hotspotAbbExonIntGr)))),
                                 fraction = c(sum(as.numeric(width(hotspotAbbExonIntGr))) / sum(as.numeric(width(hotspotAbbGr))),
                                              (sum(as.numeric(width(refGeneExonGr))) - sum(as.numeric(width(hotspotAbbExonIntGr)))) / sum(as.numeric(width(hg19Gr)))))

hotspotAbbExonOrP = poisson.test(hotspotAbbExonContMat[,1], hotspotAbbExonContMat[,2], alternative = "two.sided")
hotspotAbbExonOrR = hotspotAbbExonContMat[,1] / hotspotAbbExonContMat[,2]
hotspotAbbExonOrDf = data.frame(feature = "Exon",
                                odds_ratio = hotspotAbbExonOrR[1] / hotspotAbbExonOrR[2],
                                poisson_p = hotspotAbbExonOrP$p.value)


hotspotAbbSpliceIntGr = intersect(hotspotAbbGr, refGeneSpliceGr)
hotspotAbbSpliceContMat = matrix(c(sum(as.numeric(width(hotspotAbbSpliceIntGr))),
                                sum(as.numeric(width(refGeneSpliceGr))) - sum(as.numeric(width(hotspotAbbSpliceIntGr))),
                                sum(as.numeric(width(hotspotAbbGr))) - sum(as.numeric(width(hotspotAbbSpliceIntGr))),
                                sum(as.numeric(width(hg19Gr))) - sum(as.numeric(width(hotspotAbbGr))) - sum(as.numeric(width(refGeneSpliceGr))) + sum(as.numeric(width(hotspotAbbSpliceIntGr)))),
                              nrow = 2, dimnames = list(Hotspot = c("Hotpot", "Non-hotspot"), Feature = c("Splice", "Non-slice")))
hotspotAbbSpliceIntG = g.test(hotspotAbbSpliceContMat)

hotspotAbbSpliceIntDf = data.frame(hotspot = c("Hotspot", "Non-hotspot"),
                                 feature = c("Splice", "Splice"),
                                 size = c(sum(as.numeric(width(hotspotAbbSpliceIntGr))),
                                           sum(as.numeric(width(refGeneSpliceGr))) - sum(as.numeric(width(hotspotAbbSpliceIntGr)))),
                                 fraction = c(sum(as.numeric(width(hotspotAbbSpliceIntGr))) / sum(as.numeric(width(hotspotAbbGr))),
                                              (sum(as.numeric(width(refGeneSpliceGr))) - sum(as.numeric(width(hotspotAbbSpliceIntGr)))) / sum(as.numeric(width(hg19Gr)))))

hotspotAbbSpliceOrP = poisson.test(hotspotAbbSpliceContMat[,1], hotspotAbbSpliceContMat[,2], alternative = "two.sided")
hotspotAbbSpliceOrR = hotspotAbbSpliceContMat[,1] / hotspotAbbSpliceContMat[,2]
hotspotAbbSpliceOrDf = data.frame(feature = "Splice",
                                odds_ratio = hotspotAbbSpliceOrR[1] / hotspotAbbSpliceOrR[2],
                                poisson_p = hotspotAbbSpliceOrP$p.value)

#hotspotAbbGenicIntGr = intersect(hotspotAbbGr, refGeneGenicGr)
#hotspotAbbGenicContMat = matrix(c(sum(as.numeric(width(hotspotAbbGenicIntGr))),
                                #sum(as.numeric(width(refGeneGenicGr))) - sum(as.numeric(width(hotspotAbbGenicIntGr))),
                                #sum(as.numeric(width(hotspotAbbGr))) - sum(as.numeric(width(hotspotAbbGenicIntGr))),
                                #sum(as.numeric(width(hg19Gr))) - sum(as.numeric(width(hotspotAbbGr))) - sum(as.numeric(width(refGeneGenicGr))) + sum(as.numeric(width(hotspotAbbGenicIntGr)))),
                              #nrow = 2, dimnames = list(Hotspot = c("Hotpot", "Non-hotspot"), Feature = c("Transcript", "Non-transcript")))
#hotspotAbbGenicIntG = g.test(hotspotAbbGenicContMat)

#hotspotAbbGenicIntDf = data.frame(hotspot = c("Hotspot", "Non-hotspot"),
                                 #feature = c("Genic", "Genic"),
                                 #size = c(sum(as.numeric(width(hotspotAbbGenicIntGr))),
                                           #sum(as.numeric(width(refGeneGenicGr))) - sum(as.numeric(width(hotspotAbbGenicIntGr)))),
                                 #fraction = c(sum(as.numeric(width(hotspotAbbGenicIntGr))) / sum(as.numeric(width(hotspotAbbGr))),
                                              #(sum(as.numeric(width(refGeneGenicGr))) - sum(as.numeric(width(hotspotAbbGenicIntGr)))) / sum(as.numeric(width(hg19Gr)))))


hotspotAbbIntronIntGr = intersect(hotspotAbbGr, refGeneIntronGr)
hotspotAbbIntronContMat = matrix(c(sum(as.numeric(width(hotspotAbbIntronIntGr))),
                                sum(as.numeric(width(refGeneIntronGr))) - sum(as.numeric(width(hotspotAbbIntronIntGr))),
                                sum(as.numeric(width(hotspotAbbGr))) - sum(as.numeric(width(hotspotAbbIntronIntGr))),
                                sum(as.numeric(width(hg19Gr))) - sum(as.numeric(width(hotspotAbbGr))) - sum(as.numeric(width(refGeneIntronGr))) + sum(as.numeric(width(hotspotAbbIntronIntGr)))),
                              nrow = 2, dimnames = list(Hotspot = c("Hotpot", "Non-hotspot"), Feature = c("Intronic", "Non-intronic")))
hotspotAbbIntronIntG = g.test(hotspotAbbIntronContMat)

hotspotAbbIntronIntDf = data.frame(hotspot = c("Hotspot", "Non-hotspot"),
                                 feature = c("Intron", "Intron"),
                                 size = c(sum(as.numeric(width(hotspotAbbIntronIntGr))),
                                           sum(as.numeric(width(refGeneIntronGr))) - sum(as.numeric(width(hotspotAbbIntronIntGr)))),
                                 fraction = c(sum(as.numeric(width(hotspotAbbIntronIntGr))) / sum(as.numeric(width(hotspotAbbGr))),
                                              (sum(as.numeric(width(refGeneIntronGr))) - sum(as.numeric(width(hotspotAbbIntronIntGr)))) / sum(as.numeric(width(hg19Gr)))))

hotspotAbbIntronOrP = poisson.test(hotspotAbbIntronContMat[,1], hotspotAbbIntronContMat[,2], alternative = "two.sided")
hotspotAbbIntronOrR = hotspotAbbIntronContMat[,1] / hotspotAbbIntronContMat[,2]
hotspotAbbIntronOrDf = data.frame(feature = "Intron",
                               odds_ratio = hotspotAbbIntronOrR[1] / hotspotAbbIntronOrR[2],
                               poisson_p = hotspotAbbIntronOrP$p.value)


hotspotAbbIgrIntGr = intersect(hotspotAbbGr, refGeneIgrGr)
hotspotAbbIgrContMat = matrix(c(sum(as.numeric(width(hotspotAbbIgrIntGr))),
                                sum(as.numeric(width(refGeneIgrGr))) - sum(as.numeric(width(hotspotAbbIgrIntGr))),
                                sum(as.numeric(width(hotspotAbbGr))) - sum(as.numeric(width(hotspotAbbIgrIntGr))),
                                sum(as.numeric(width(hg19Gr))) - sum(as.numeric(width(hotspotAbbGr))) - sum(as.numeric(width(refGeneIgrGr))) + sum(as.numeric(width(hotspotAbbIgrIntGr)))),
                              nrow = 2, dimnames = list(Hotspot = c("Hotpot", "Non-hotspot"), Feature = c("Intergenic", "Non-integenic")))
hotspotAbbIgrIntG = g.test(hotspotAbbIgrContMat)

hotspotAbbIgrIntDf = data.frame(hotspot = c("Hotspot", "Non-hotspot"),
                                 feature = c("Intergenic", "Intergenic"),
                                 size = c(sum(as.numeric(width(hotspotAbbIgrIntGr))),
                                           sum(as.numeric(width(refGeneIgrGr))) - sum(as.numeric(width(hotspotAbbIgrIntGr)))),
                                 fraction = c(sum(as.numeric(width(hotspotAbbIgrIntGr))) / sum(as.numeric(width(hotspotAbbGr))),
                                              (sum(as.numeric(width(refGeneIgrGr))) - sum(as.numeric(width(hotspotAbbIgrIntGr)))) / sum(as.numeric(width(hg19Gr)))))

hotspotAbbIgrOrP = poisson.test(hotspotAbbIgrContMat[,1], hotspotAbbIgrContMat[,2], alternative = "two.sided")
hotspotAbbIgrOrR = hotspotAbbIgrContMat[,1] / hotspotAbbIgrContMat[,2]
hotspotAbbIgrOrDf = data.frame(feature = "Intergenic",
                               odds_ratio = hotspotAbbIgrOrR[1] / hotspotAbbIgrOrR[2],
                               poisson_p = hotspotAbbIgrOrP$p.value)




#
hotspotAbbIntDf = rbind(hotspotAbb3UtrIntDf, hotspotAbb5UtrIntDf, hotspotAbbExonIntDf, hotspotAbbIntronIntDf, hotspotAbbIgrIntDf)

p = ggplot(data = hotspotAbbIntDf, aes(x=hotspot, y=fraction)) +
    geom_bar(stat="identity") +
    theme(plot.title   = element_text(size = baseFontSize, face="bold"),
          axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90), axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
          axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 45, hjust=1),
          axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
          strip.text.x = element_text(size = baseFontSize, face="bold"),
          strip.text.y = element_text(size = baseFontSize, face="bold", angle = 0),
          legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
          legend.text  = element_text(size = baseFontSize, family="sans"),
          legend.direction = "horizontal",
          legend.position = "none") +
    scale_x_discrete(name="") +
    scale_y_continuous("Percentage\n", label = percent) +
    facet_wrap(~ feature, scales = "free_y", nrow = 1)

hotspotAbbIntFile = file.path(baseDir, "figure", "Genomic_enrichment_of_hotspot100_fdr0.01.pdf")
ggsave(filename = hotspotAbbIntFile, plot = p, width = 18, height = 8, unit = "cm")

#
hotspotAbbOrDf = rbind(hotspotAbb3UtrOrDf, hotspotAbb5UtrOrDf, hotspotAbbExonOrDf, hotspotAbbIntronOrDf, hotspotAbbIgrOrDf)

p = ggplot(data = hotspotAbbOrDf, aes(x=feature, y=log(odds_ratio))) +
    geom_bar(stat="identity") +
    theme(plot.title   = element_text(size = baseFontSize, face="bold"),
          axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90), axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
          axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 45, hjust=1),
          axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
          strip.text.x = element_text(size = baseFontSize, face="bold"),
          strip.text.y = element_text(size = baseFontSize, face="bold", angle = 0),
          legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
          legend.text  = element_text(size = baseFontSize, family="sans"),
          legend.direction = "horizontal",
          legend.position = "none") +
    scale_x_discrete(name="") +
    scale_y_continuous("Log2 odds ratio\n")

hotspotAbbOrFile = file.path(baseDir, "figure", "Genomic_enrichment_of_hotspot100_fdr0.01.odds.pdf")
ggsave(filename = hotspotAbbOrFile, plot = p, width = 8, height = 8, unit = "cm")


# Epigenomic features
chromHmmFiles = Sys.glob(file.path(baseDir, "roadmap/chromhmm/*_coreMarks_mnemonics.bed.gz"))
chromHmmAllDf = data.frame()
for (chromHmmFile in chromHmmFiles) {
    cat(sprintf("Reading %s ...\n", chromHmmFile))
    tissueType = strsplit(basename(chromHmmFile), "_", fixed = T)[[1]][1]
    chromHmmDf = read.delim(gzfile(chromHmmFile), header = F, as.is = T)
    colnames(chromHmmDf) = c("space", "start", "end", "state")
    chromHmmDf$type = tissueType
    chromHmmAllDf = rbind(chromHmmAllDf, chromHmmDf)
}

chromHmmAllDf$space = gsub("chr", "", chromHmmAllDf$space)
chromHmmAllDf = chromHmmAllDf[chromHmmAllDf$space %in% chrs,]
chromHmmAllDf$space = factor(chromHmmAllDf$space, levels = chrs)

hotspotAbbStateIntAllDf = data.frame()
hotspotAbbStateOrAllDf = data.frame()

for (s in mixedsort(unique(chromHmmAllDf$state))) {
    print(s)
    chromHmmStateAllDf = subset(chromHmmAllDf, state == s)
    chromHmmStateAllGr = reduce(as(as(chromHmmStateAllDf, "RangedData"), "GRanges"))
    hotspotAbbStateIntGr = intersect(hotspotAbbGr, chromHmmStateAllGr)
    hotspotAbbStateContMat = matrix(c(sum(as.numeric(width(hotspotAbbStateIntGr))),
                                    sum(as.numeric(width(chromHmmStateAllGr))) - sum(as.numeric(width(hotspotAbbStateIntGr))),
                                    sum(as.numeric(width(hotspotAbbGr))) - sum(as.numeric(width(hotspotAbbStateIntGr))),
                                    sum(as.numeric(width(hg19Gr))) - sum(as.numeric(width(hotspotAbbGr))) - sum(as.numeric(width(chromHmmStateAllGr))) + sum(as.numeric(width(hotspotAbbStateIntGr)))),
                                  nrow = 2, dimnames = list(Hotspot = c("Hotpot", "Non-hotspot"), Feature = c(s, sprintf("Non-%s", s))))
    hotspotAbbStateIntG = g.test(hotspotAbbStateContMat)
    #hotspotAbbStateIntC = chisq.test(hotspotAbbStateContMat)

    hotspotAbbStateIntDf = data.frame(hotspot = c("Hotspot", "Non-hotspot"),
                                     feature = c(s, s),
                                     size = c(sum(as.numeric(width(hotspotAbbStateIntGr))),
                                               sum(as.numeric(width(chromHmmStateAllGr))) - sum(as.numeric(width(hotspotAbbStateIntGr)))),
                                     fraction = c(sum(as.numeric(width(hotspotAbbStateIntGr))) / sum(as.numeric(width(hotspotAbbGr))),
                                                  (sum(as.numeric(width(chromHmmStateAllGr))) - sum(as.numeric(width(hotspotAbbStateIntGr)))) / sum(as.numeric(width(hg19Gr)))))
    hotspotAbbStateIntAllDf = rbind(hotspotAbbStateIntAllDf, hotspotAbbStateIntDf)

    #hotspotAbbStateOrP = poisson.test(hotspotAbbStateContMat[,1], hotspotAbbStateContMat[,2], alternative = "two.sided")
    hotspotAbbStateOrR = hotspotAbbStateContMat[,1] / hotspotAbbStateContMat[,2]
    hotspotAbbStateOrDf = data.frame(feature = s,
                                     odds_ratio = hotspotAbbStateOrR[1] / hotspotAbbStateOrR[2],
                                     gtest_p = hotspotAbbStateIntG$p.value)

    hotspotAbbStateOrAllDf = rbind(hotspotAbbStateOrAllDf, hotspotAbbStateOrDf)
}

#
p = ggplot(data = hotspotAbbStateIntAllDf, aes(x=hotspot, y=fraction, color = feature, fill = feature)) +
    geom_bar(stat="identity") +
    theme(plot.title   = element_text(size = baseFontSize, face="bold"),
          axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
          axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
          axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 45, hjust=1),
          axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
          strip.text.x = element_text(size = baseFontSize, face="bold"),
          strip.text.y = element_text(size = baseFontSize, face="bold", angle = 0),
          legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
          legend.text  = element_text(size = baseFontSize, family="sans"),
          legend.direction = "horizontal",
          legend.position = "none") +
    scale_x_discrete(name="") +
    scale_y_continuous("Percentage\n", label = percent) +
    scale_color_manual("", values = roadmapChromBorderColors) +
    scale_fill_manual("", values = roadmapChromFillColors) +
    facet_wrap(~ feature, scales = "free_y", ncol = 5)

hotspotAbbStateIntAllFile = file.path(baseDir, "figure", "Epignomic_enrichment_of_hotspot100_fdr0.01.pdf")
ggsave(filename = hotspotAbbStateIntAllFile, plot = p, width = 28, height = 28, unit = "cm")

#
p = ggplot(data = hotspotAbbStateOrAllDf, aes(x=feature, y=log2(odds_ratio), color = feature, fill = feature)) +
    geom_bar(stat="identity") +
    theme(plot.title   = element_text(size = baseFontSize, face="bold"),
          axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
          axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
          axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 45, hjust=1),
          axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
          strip.text.x = element_text(size = baseFontSize, face="bold"),
          strip.text.y = element_text(size = baseFontSize, face="bold", angle = 0),
          legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
          legend.text  = element_text(size = baseFontSize, family="sans"),
          legend.direction = "horizontal",
          legend.position = "none") +
    scale_x_discrete(name="") +
    scale_y_continuous("Log2 odds ratio\n") +
    scale_color_manual("", values = roadmapChromBorderColors) +
    scale_fill_manual("", values = roadmapChromFillColors)

hotspotAbbStateOrAllFile = file.path(baseDir, "figure", "Epignomic_enrichment_of_hotspot100_fdr0.01.odds.pdf")
ggsave(filename = hotspotAbbStateOrAllFile, plot = p, width = 14, height = 14, unit = "cm")



##
##
##
hotspotGenomicDistDf = data.frame()
hotspotEpigenomicDistDf = data.frame()
for (chr in chrs) {
    cat(sprintf("Processing chr%s ...\n", chr))
    hotspotChrLabelDf = subset(hotspotSigAssAnnotDf, space == chr)
    for (i in 1:nrow(hotspotChrLabelDf)) {
        hotspot = hotspotChrLabelDf[i,]
        segStates = unique(strsplit(hotspot$state, ",")[[1]])
        segStates = segStates[segStates != ""]
        if (length(segStates) == 0) segStates = NA
        hotspotEpigenomicDistDf = rbind(hotspotEpigenomicDistDf,
                                        data.frame(Location = hotspot$Location,
                                                   all_adj_mcnt_p_value_d = hotspot$all_adj_mcnt_p_value_d,
                                                   state = segStates,
                                                   stringsAsFactors = FALSE))
        hotspotSigAnnotSubDf = subset(hotspotSigAnnotDf, Location == hotspot$Location)
        splice_site_cnt = sum(grepl("splice", hotspot$types))
        exon_cnt = sum(grepl("exon", hotspot$types))
        five_prime_utr_cnt = sum(grepl("five", hotspot$types))
        three_prime_utr_cnt = sum(grepl("three", hotspot$types))
        intron_cnt = sum(grepl("intron", hotspot$types))
        tot_region_cnt = splice_site_cnt + exon_cnt + five_prime_utr_cnt + three_prime_utr_cnt + intron_cnt
        if (splice_site_cnt > 0) hotspotGenomicDistDf = rbind(hotspotGenomicDistDf,
                                                              data.frame(Location = hotspot$Location,
                                                                         all_adj_mcnt_p_value_d = hotspot$all_adj_mcnt_p_value_d,
                                                                         region = "splice_site"))
        if (exon_cnt > 0) hotspotGenomicDistDf = rbind(hotspotGenomicDistDf, 
                                                       data.frame(Location = hotspot$Location, 
                                                                  all_adj_mcnt_p_value_d = hotspot$all_adj_mcnt_p_value_d,
                                                                  region = "exon"))
        if (five_prime_utr_cnt > 0) hotspotGenomicDistDf = rbind(hotspotGenomicDistDf, 
                                                                 data.frame(Location = hotspot$Location, 
                                                                            all_adj_mcnt_p_value_d = hotspot$all_adj_mcnt_p_value_d, 
                                                                            region = "five_prime_utr"))
        if (three_prime_utr_cnt > 0) hotspotGenomicDistDf = rbind(hotspotGenomicDistDf, 
                                                                  data.frame(Location = hotspot$Location, 
                                                                             all_adj_mcnt_p_value_d = hotspot$all_adj_mcnt_p_value_d, 
                                                                             region = "three_prime_utr"))
        if (intron_cnt > 0) hotspotGenomicDistDf = rbind(hotspotGenomicDistDf, 
                                                         data.frame(Location = hotspot$Location, 
                                                                    all_adj_mcnt_p_value_d = hotspot$all_adj_mcnt_p_value_d, 
                                                                    region = "intron"))
        if (tot_region_cnt == 0) hotspotGenomicDistDf = rbind(hotspotGenomicDistDf, 
                                                              data.frame(Location = hotspot$Location, 
                                                                         all_adj_mcnt_p_value_d = hotspot$all_adj_mcnt_p_value_d, 
                                                                         region = "igr"))
    }
}

hotspotGenomicDistDf = hotspotGenomicDistDf[!grepl("^Y", hotspotGenomicDistDf$Location, perl = T),]
hotspotEpigenomicDistDf = hotspotEpigenomicDistDf[!grepl("^Y", hotspotEpigenomicDistDf$Location, perl = T),]
hotspotGenomicDistGrpDf = sqldf('SELECT region, COUNT(*) AS cnt FROM hotspotGenomicDistDf GROUP BY region')
hotspotGenomicDistGrpDf$freq = hotspotGenomicDistGrpDf$cnt / length(unique(hotspotGenomicDistDf$Location))
hotspotEpigenomicDistGrpDf = sqldf('SELECT state, COUNT(*) AS cnt FROM hotspotEpigenomicDistDf GROUP BY state')
hotspotEpigenomicDistGrpDf$freq = hotspotEpigenomicDistGrpDf$cnt / length(unique(hotspotEpigenomicDistDf$Location))
hotspotEpigenomicDistGrpDf[1, "state"] = "Unk"


##
## Plot genomic and epigenomic distributions of hotspots
##
baseFontSize = 11
hotspotGenomicDistGrpDf = hotspotGenomicDistGrpDf[order(-hotspotGenomicDistGrpDf$cnt),]
p = ggplot(data = hotspotGenomicDistGrpDf, aes(x=factor(region, levels = hotspotGenomicDistGrpDf$region), y=freq)) +
    geom_bar(stat="identity") +
    theme(plot.title   = element_text(size = baseFontSize, face="bold"),
          axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
          axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
          axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 45, hjust=1),
          axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
          strip.text.x = element_text(size = baseFontSize, face="bold"),
          strip.text.y = element_text(size = baseFontSize, face="bold", angle = 0),
          legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
          legend.text  = element_text(size = baseFontSize, family="sans"),
          legend.direction = "horizontal",
          legend.position = "none") +
    scale_x_discrete(name="\nGenomic region", labels = c("intron" = "Intron", "igr" = "Intergenic region",
                                                         "splice_site" = "Splice region", "five_prime_utr" = "5' UTR",
                                                         "three_prime_utr" = "3' UTR", "exon" = "Exon")) +
    #scale_x_discrete(name="\nGenomic region") +
    scale_y_continuous("Percentage\n", label = percent)
hotspotGenomicDistGrpFile = file.path(baseDir, "figure", "Genomic_distribution_of_hotspot100_fdr0.01.pdf")
ggsave(filename = hotspotGenomicDistGrpFile, plot = p, width = 8, height = 8, unit = "cm")

baseFontSize = 11
hotspotEpigenomicDistGrpDf = hotspotEpigenomicDistGrpDf[order(-hotspotEpigenomicDistGrpDf$cnt),]
p = ggplot(data = hotspotEpigenomicDistGrpDf, aes(x=factor(state, levels = hotspotEpigenomicDistGrpDf$state), y=freq, fill = state, color = state)) +
    geom_bar(stat="identity") +
    theme(plot.title   = element_text(size = baseFontSize, face="bold"),
          axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
          axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
          axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 45, hjust=1),
          axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
          strip.text.x = element_text(size = baseFontSize, face="bold"),
          strip.text.y = element_text(size = baseFontSize, face="bold", angle = 0),
          legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
          legend.text  = element_text(size = baseFontSize, family="sans"),
          legend.direction = "horizontal",
          legend.position = "none") +
    scale_x_discrete(name="\nChromatin state", labels = roadmapChromStateLabels) +
    scale_y_continuous("Percentage\n", label = percent) +
    scale_color_manual("", values = roadmapChromBorderColors) +
    scale_fill_manual("", values = roadmapChromFillColors)
hotspotEpigenomicDistGrpFile = file.path(baseDir, "figure", "Epigenomic_distribution_of_hotspot100_fdr0.01.pdf")
ggsave(filename = hotspotEpigenomicDistGrpFile, plot = p, width = 14, height = 12, unit = "cm")


##
## Correlation analysis between mutation rates and epigenomic features
##

hotspotSigAnnotGr = reduce(sort(as(as(hotspotSigAnnotDf, "RangedData"), "GRanges")))

cancerTypes = unique(wgsSnvAllDf$cancer)

cancerSnvDf = subset(wgsSnvAllDf, cancer == "LICA")
cancerSnvDf$start = cancerSnvDf$pos
cancerSnvDf$end = cancerSnvDf$pos
cancerSnvDf$space = factor(cancerSnvDf$space, levels=chrs)
cancerSnvGr = sort(as(as(cancerSnvDf, "RangedData"), "GRanges"))

cancerHotspotSnvGr = cancerSnvGr[cancerSnvGr %over% hotspotSigAnnotGr]
cancerNonHotspotSnvGr = cancerSnvGr[cancerSnvGr %outside% hotspotSigAnnotGr]

# Read Roadmap EID mapping table
eidMapFile = file.path(baseDir, "roadmap", "Roadmap_EID_Mapping_Table.txt")
eidMapDf = read.delim(eidMapFile, header = F, as.is = T)
colnames(eidMapDf) = c("eid", "group", "color", "desc")

#brainEids = subset(eidMapDf, group == "Breast")$eid
selEids = c("E066", "E118")

# Read Roadmap peak bed files
totPeakLst <- foreach(eid=selEids) %dopar% {
    print(eid)
    peakEidBedFiles = Sys.glob(file.path(baseDir, "roadmap", "peaks", "narrowPeak", sprintf("%s-*.narrowPeak.gz", eid)))
    totEidPeakLst = list()
    for (peakEidBedFile in peakEidBedFiles) {
        peakType = strsplit(strsplit(basename(peakEidBedFile), "-")[[1]][2], ".", fixed = T)[[1]][1]
        print(peakType)
        peakBedDf = read.delim(gzfile(peakEidBedFile), header = F, as.is = T)[, c(1:3)]
        colnames(peakBedDf) = c("space", "start", "end")
        peakBedDf$type = peakType
        peakBedDf$space = gsub("chr", "", peakBedDf$space)
        peakBedDf = subset(peakBedDf, space %in% chrs)
        peakBedDf$space = factor(peakBedDf$space, levels = chrs)
        peakBedDf = with(peakBedDf, peakBedDf[order(space, start, end),])
        peakBedGr = as(as(peakBedDf, "RangedData"), "GRanges")
        totEidPeakLst[[peakType]] = peakBedGr
    }
    totEidPeakLst
}
names(totPeakLst) = selEids

# Generate tiles and count SNVs in single neurons and TCGA GBM WGS
seqlens = hg19Df$size
names(seqlens) = hg19Df$chrom
tiles = tileGenome(seqlengths=seqlens, tilewidth=10^7, cut.last.tile.in.chrom=TRUE)
tiles$hotspot_snv_cnt = countOverlaps(tiles, cancerHotspotSnvGr)
tiles$non_hotspot_snv_cnt = countOverlaps(tiles, cancerNonHotspotSnvGr)

mutPeakCorrDf = data.frame()
for (j in 1:length(totPeakLst)) {
    eid = names(totPeakLst)[j]
    totEidPeakLst = totPeakLst[[j]]
    print(eid)
    mutPeakCorrEidDf <- foreach(k=1:length(totEidPeakLst), .combine = rbind) %dopar% {
        peakType = names(totEidPeakLst)[k]
        peakGr = totEidPeakLst[[k]]
        peakColName = sprintf("%s_%s", eid, peakType)
        print(peakColName)
        values(tiles)[, peakColName] = 0
        for (l in 1:length(tiles)) {
            tile = tiles[l]
            olMutPeak = findOverlaps(tile, peakGr)
            values(tiles)[l, peakColName] = sum(width(pintersect(tile[queryHits(olMutPeak)], peakGr[subjectHits(olMutPeak)])))
        }
        scorr_hotspot = cor.test(values(tiles)[, peakColName], tiles$hotspot_snv_cnt, method = "spearman")
        pcorr_hotspot = cor.test(values(tiles)[, peakColName], tiles$hotspot_snv_cnt, method = "pearson")
        scorr_non_hotspot = cor.test(values(tiles)[, peakColName], tiles$non_hotspot_snv_cnt, method = "spearman")
        pcorr_non_hotspot = cor.test(values(tiles)[, peakColName], tiles$non_hotspot_snv_cnt, method = "pearson")
        rbind(
              data.frame(eid = eid, peak = peakType,
                         sample_type = "Hotspot SNVs",
                         pcorr = as.vector(pcorr_hotspot$estimate), pp = pcorr_hotspot$p.value,
                         scorr = as.vector(scorr_hotspot$estimate), sp = scorr_hotspot$p.value),
              data.frame(eid = eid, peak = peakType,
                         sample_type = "Non-hotspot SNVs",
                         pcorr = as.vector(pcorr_non_hotspot$estimate), pp = pcorr_non_hotspot$p.value,
                         scorr = as.vector(scorr_non_hotspot$estimate), sp = scorr_non_hotspot$p.value))
    }
    mutPeakCorrDf = rbind(mutPeakCorrDf, mutPeakCorrEidDf)
}

mutPeakCorrDf = merge(mutPeakCorrDf, eidMapDf, by = c("eid"), all.x = T)

# Plot correlation
for (eid in selEids) {
    print(eid)
    mutPeakCorrEidDf = mutPeakCorrDf[mutPeakCorrDf$eid == eid,]
    mutPeakCorrEidDf = with(mutPeakCorrEidDf, mutPeakCorrEidDf[order(pcorr, decreasing = T),])
    mutPeakCorrEidSngDf = mutPeakCorrEidDf[mutPeakCorrEidDf$sample_type == "Hotspot SNVs",]
    baseFontSize = 11
    p = ggplot(data = mutPeakCorrEidDf, aes(x=factor(peak, levels = mutPeakCorrEidSngDf$peak), y=pcorr, fill = sample_type)) +
        geom_bar(stat="identity", position="dodge", width = 0.5) +
        ggtitle(sprintf("%s (%s)\n", mutPeakCorrEidDf$desc[1], eid)) +
        theme(plot.title   = element_text(size = baseFontSize + 3, face="bold"),
              plot.margin  = unit(c(1, 1, 1, 2), "cm"),
              axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
              axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
              axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 45, hjust=1),
              axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
              strip.text.x = element_text(size = baseFontSize, face="bold"),
              strip.text.y = element_text(size = baseFontSize, face="bold", angle = 0),
              legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
              legend.text  = element_text(size = baseFontSize, family="sans"),
              legend.direction = "vertical",
              legend.position = "right") +
        scale_x_discrete(name="") +
        scale_y_continuous("Pearson correlation\n") +
        scale_fill_manual(values = c("#A40000", "#005C95"), guide = guide_legend(title = NULL, ncol = 1))
        plotFile = file.path(baseDir, "figure", sprintf("%s-Correlation_between_BRCA_mutation_rates_and_Roadmap_peaks.pdf", eid))
        ggsave(filename = plotFile, plot = p, width=18, height=15, units = "cm")
}


##
## Hotspot distribution across cancer types
##
sampleIdsByCancerLst = list()
for (c in unique(wgsSnvAllGrpBySidDf$cancer)) {
    print(c)
    wgsSnvAllGrpBySidCancerDf = subset(wgsSnvAllGrpBySidDf, cancer == c)
    sampleIdsByCancerLst[[c]] = wgsSnvAllGrpBySidCancerDf$sid
}

for (i in 1:nrow(hotspotSigAnnotGrpDf)) {
    print(100 * i/nrow(hotspotSigAnnotGrpDf))
    sampleIds = strsplit(hotspotSigAnnotGrpDf[i, "sids"], ",")[[1]]
    for (c in unique(wgsSnvAllGrpBySidDf$cancer)) {
        hotspotSigAnnotGrpDf[i, sprintf("%s_CNT", c)] = sum(sampleIds %in% sampleIdsByCancerLst[[c]])
        hotspotSigAnnotGrpDf[i, sprintf("%s_PCT", c)] = 100 * sum(sampleIds %in% sampleIdsByCancerLst[[c]]) / length(sampleIdsByCancerLst[[c]])
        hotspotSigAnnotGrpDf[i, sprintf("%s_BIN", c)] = sum(sampleIds %in% sampleIdsByCancerLst[[c]]) > 0
    }
}

hotspotSigAnnotGrpDf = hotspotSigAnnotGrpDf[order(hotspotSigAnnotGrpDf$all_adj_scnt_p_value),]
hotspotSigAnnotGrpSubDf = subset(hotspotSigAnnotGrpDf, all_adj_mcnt_p_value < 0.001 & all_adj_scnt_p_value < 0.001)
hotspotSigAnnotGrpSubDf = hotspotSigAnnotGrpSubDf[, -which(colnames(hotspotSigAnnotGrpSubDf) %in% c(grep("_CNT", colnames(hotspotSigAnnotGrpDf), value = T), grep("_BIN", colnames(hotspotSigAnnotGrpDf), value = T)))]

nrow(hotspotSigAnnotGrpSubDf)
hotspotSigAnnotGrpSubFile = file.path(tableDir, "Distirution_of_hotspots_in_cancer_types.txt")
write.table(hotspotSigAnnotGrpSubDf, hotspotSigAnnotGrpSubFile, row.names = F, col.names = T, sep = "\t", quote = F)

hotspotSigAnnotGrpPctDf = hotspotSigAnnotGrpDf[, c(grep("_PCT", colnames(hotspotSigAnnotGrpDf)))]
hotspotSigAnnotGrpMat = log10(as.matrix(head(hotspotSigAnnotGrpDf[, c(grep("_PCT", colnames(hotspotSigAnnotGrpDf)))], 50)))
hotspotSigAnnotGrpMat = hotspotSigAnnotGrpMat[rev(rownames(hotspotSigAnnotGrpMat)),]

colfunc = colorRampPalette(c("white", "red"))
image(t(hotspotSigAnnotGrpMat), col = colfunc(100), axes=FALSE, xlab="", ylab="")
axis(3, at = 1:ncol(hotspotSigAnnotGrpMat), labels=colnames(hotspotSigAnnotGrpMat), tick=T)
axis(2, at = 1:nrow(hotspotSigAnnotGrpMat), labels=rownames(hotspotSigAnnotGrpMat), tick=T)



##
## Examples of hotspots from each genomic/epigenomic category
##

## Exonic (Active Transcript)
hotspotExonicDf = subset(hotspotSigAssAnnotDf, grepl("exon", types) & !grepl("utr", types))
nrow(hotspotExonicDf)
names(hotspotExonicDf)
head(hotspotExonicDf[, c(1, 5, 6, 8:11, 14)], 20)
tail(hotspotExonicDf[, c(1, 5, 6, 8:11, 14)], 10)

#write.table(head(hotspotSigAssAnnotDf[, c("Location", "mcnt", "scnt", "all_adj_mcnt_p_value", "all_adj_scnt_p_value", "overlapping_genes", "cis_genes", "states")], 1000),
write.table(head(hotspotSigAssAnnotDf, 1000),
            file = "hotspot1000.txt", row.names = F, col.names = T, sep = "\t", quote = F)

#loc = "X:18122607-18122755"
#hotspotSubDf = subset(hotspotExonicDf, Location == loc)
nrow(hotspotSigAssAnnotDf)

gene = "CTNNB1"
testHotDf = subset(hotspotSigAssAnnotDf, grepl(gene, cis_genes) | grepl(gene, overlapping_genes))
head(testHotDf[, c("Location", "mcnt", "scnt", "all_adj_mcnt_p_value", "all_adj_scnt_p_value", "overlapping_genes", "cis_genes", "states")], 20)

## Promoter (TSS)

##
## Trans-interactions bewteen distal elements and promoters
##
tssDf = subset(hotspotSigAssAnnotDf, cis_genes != "" & grepl("Tss", states))
nrow(tssDf)
head(tssDf[, c("Location", "mcnt", "scnt", "all_adj_mcnt_p_value", "all_adj_scnt_p_value", "overlapping_genes", "cis_genes", "states")], 20)

enhDf = subset(hotspotSigAssAnnotDf, cis_genes != "" & grepl("Enh", states) & !grepl("Tss", states))
nrow(enhDf)
head(enhDf[, c("Location", "mcnt", "scnt", "all_adj_mcnt_p_value", "all_adj_scnt_p_value", "overlapping_genes", "cis_genes", "states")], 20)

ctcfDf = subset(hotspotSigAssAnnotDf, cis_genes != "" & grepl("ZNF", states) & !grepl("Tss", states))
nrow(ctcfDf)
head(ctcfDf[, c("Location", "mcnt", "scnt", "all_adj_mcnt_p_value", "all_adj_scnt_p_value", "overlapping_genes", "cis_genes", "states")], 20)


##
## Cancer distribution analysis
##
#cancerToColorDf = data.frame(cancer = c("BRCA", "COAD", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "OV", "READ", "SKCM", "STAD", "THCA", "UCEC"),
                             #color = c("#000000", "#0B2E82", "#006A3A", "#663E25", "#6F0042", "#913696", "#CC9DCB", "#B48541", "#CA012D", "#FF3600", "#B26345", "#01874A", "#DF54A7", "#FCBA64", "#FFC211", "#4FB0D1"))

#vcfFiles = Sys.glob(file.path(vcfDir, "cancer/*/*.stage4.vcf.gz"))
#sampleIdsByCancerLst = list()
#sampleIdToCancerDf = data.frame()
#for (vcfFile in vcfFiles) {
    #cancerType = basename(dirname(vcfFile))
    #sampleIds = samples(scanVcfHeader(vcfFile))
    #sampleIdsByCancerLst[[cancerType]] = sampleIds
    #for (sampleId in sampleIds) {
        #sampleIdToCancerDf = rbind(sampleIdToCancerDf, data.frame(id = sampleId, cancer = cancerType))
    #}
#}

#totSampleIds = as.vector(unlist(sampleIdsByCancerLst))
#totAbbSampleIds = as.vector(sapply(totSampleIds, function(x) paste(strsplit(x, "-")[[1]][2:3], collapse = "-")))
#sampleIdToCancerDf$sid = sapply(as.character(sampleIdToCancerDf$id), function(x) paste(strsplit(x, "-")[[1]][2:3], collapse = "-"))


## Create color range file for ITOL
# http://itol.embl.de/help/help.shtml
sampleIdToCancerDf = merge(sampleIdToCancerDf, cancerToColorDf, by = "cancer")
sampleIdToCancerDf$type = "range"
sampleIdToCancerFile = file.path(baseDir, "table", "NJ_tree-hotspot100-fdr0.01.colors.txt")
write.table(sampleIdToCancerDf[, c("sid", "type", "color", "cancer")], sampleIdToCancerFile, row.names = F, col.names = F, sep = "\t", quote = F)


## Read annotated hotspot files
hotspotChrAnnotFiles = mixedsort(Sys.glob(file.path(hotspotDir, sprintf("cancer*.chr*.hotspot%d.fdr0.01.annotated.txt", distCut))))
hotspotSigAnnotDf = data.frame()
hotspotSigAnnotGrpDf = data.frame()
for (hotspotChrAnnotFile in hotspotChrAnnotFiles) {
    cat(sprintf("Reading %s ...\n", hotspotChrAnnotFile))
    hotspotChrAnnotDf = read.delim(hotspotChrAnnotFile, header = T, as.is = T)
    hotspotSigAnnotDf = rbind(hotspotSigAnnotDf, hotspotChrAnnotDf)
    hotspotChrAnnotDf$allele = NULL
    hotspotChrAnnotDf$strand = NULL
    hotspotChrAnnotGrpDf = sqldf('SELECT Location, all_adj_mcnt_p_value, all_adj_scnt_p_value, sids
                                 FROM hotspotChrAnnotDf 
                                 GROUP BY Location ORDER BY all_adj_mcnt_p_value ASC')
    hotspotSigAnnotGrpDf = rbind(hotspotSigAnnotGrpDf, hotspotChrAnnotGrpDf)
}


## Filter very significant hotspots (FDR < 0.01)
hotspotSigAnnotDf = hotspotSigAnnotDf[!grepl("^Y", hotspotSigAnnotDf$Location, perl = T),]
hotspotSigAnnotGrpDf = hotspotSigAnnotGrpDf[!grepl("^Y", hotspotSigAnnotGrpDf$Location, perl = T),]
#hotspotAnnotVerySigDf = subset(hotspotSigAnnotDf, all_adj_mcnt_p_value < 0.01)
#hotspotAnnotVerySigGrpDf = subset(hotspotSigAnnotGrpDf, all_adj_mcnt_p_value < 0.01)

## Build phylogenetic tree from distances based on hotspot profiles across samples
sampleHotspotMutProfileDf <- foreach(sampleId=totAbbSampleIds, .combine = rbind) %dopar% {
    print(sampleId)
    data.frame(t(sapply(1:nrow(hotspotSigAnnotGrpDf), 
                        function(x) as.integer(grepl(sampleId, hotspotSigAnnotGrpDf[x,]$sids)))), 
               row.names = c(sampleId))
}

sampleHotspotMutProfileDist = vegdist(sampleHotspotMutProfileDf, method = "jaccard")
#sampleHotspotMutProfileNjTree = nj(as.dist(sampleHotspotMutProfileDist))
#sampleHotspotMutProfileNjTree = nj(sampleHotspotMutProfileDist)
sampleHotspotMutProfileNjTree = hclust(sampleHotspotMutProfileDist, method = "single")
sampleHotspotMutProfileNjTreeFile = file.path(baseDir, "table", "NJ_tree-hotspot100-fdr0.01.jaccard.single.newick.txt")
write.tree(as.phylo(sampleHotspotMutProfileNjTree), file=sampleHotspotMutProfileNjTreeFile)


##
## Add annotation for known cancer genes
##
allOncoFile = file.path(baseDir, "allonco/allonco_20130923.tsv")
allOncoDf = read.delim(allOncoFile, header = T, as.is = T)
hotspotSigAnnotDf$allOnco = ifelse(hotspotSigAnnotDf$SYMBOL %in% allOncoDf$symbol, "Y", "N")

intogenDriverSignalFile = file.path(baseDir, "intogen/Mutational_drivers_signals.tsv")
intogenDriverSignalDf = read.delim(intogenDriverSignalFile, header = T, as.is = T, comment.char = "#")
hotspotSigAnnotDf$intogenDriverSignals = ifelse(hotspotSigAnnotDf$SYMBOL %in% intogenDriverSignalDf$geneHGNCsymbol, "Y", "N")

intogenDriverFile = file.path(baseDir, "intogen/Drivers_type_role.tsv")
intogenDriverDf = read.delim(intogenDriverFile, header = T, as.is = T, comment.char = "#")
hotspotSigAnnotDf$intogenDrivers = ifelse(hotspotSigAnnotDf$SYMBOL %in% intogenDriverDf$geneHGNCsymbol, "Y", "N")

cgcFile = file.path(baseDir, "cosmic/cancer_gene_census.csv")
cgcDf = read.csv(cgcFile, header = T, as.is = T)
hotspotSigAnnotDf$cosmic = ifelse(hotspotSigAnnotDf$SYMBOL %in% cgcDf$Gene.Symbol, "Y", "N")

tcgaSmgFile = file.path(baseDir, "tcga-pancan/hcd.txt")
tcgaSmgDf = read.delim(tcgaSmgFile, header = T, as.is = T)
hotspotSigAnnotDf$TCGA_PanCan_SMG = ifelse(hotspotSigAnnotDf$SYMBOL %in% tcgaSmgDf$Gene.Symbol, "Y", "N")

