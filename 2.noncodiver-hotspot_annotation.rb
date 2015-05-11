#!/usr/bin/env ruby

require 'pathname'

distCut = 100
baseDir = Pathname.new "/n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot"
hotspotChrVepSigFiles = Pathname.glob(baseDir + "hotspot" + "Combined*hotspot#{distCut}.fdr0.05.vep_in.txt")
hotspotChrVepSigFiles.each do |hotspotChrVepSigFile|
  hotspotChrVepSigOutFile = hotspotChrVepSigFile.dirname + hotspotChrVepSigFile.basename(".txt").sub_ext(".vep_out.txt")
  hotspotChrVepSigLsfOutFile = hotspotChrVepSigOutFile.sub_ext(".txt.lsfout")
  next if hotspotChrVepSigLsfOutFile.exist?
  hotspotChrVepSigLsfShFile = hotspotChrVepSigFile.sub_ext(".sh")
  hotspotChrVepSigLsfShFile.open('w') do |file|
    #--plugin UpDownDistance,1000000,1000000 \\
    content =<<-CONTENT
#!/bin/sh
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
    --fasta /home/sl279/.vep/homo_sapiens/78_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/ucsc/encode/distalDhsToPromoterDhs.bed.gz,distalDhsToPromoterDhs,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/ucsc/encode/dhsToGeneExpression.bed.gz,dhsToGeneExpression,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/fantom5/fantom5EnhancerTssAssociations.bed.gz,fantom5EnhancerTssAssociations,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/dbsuper/all/allDbSuperEnhancerGeneAssociations.bed.gz,allDbSuperEnhancerGeneAssociations,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/4dgenome/4dGenomeHomoSapiens.bed.gz,fourDGenomeHomoSapiens,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/insitu-hic/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist.annotated.bed.gz,GSE63525_GM12878_HiCCUPS_looplist,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/insitu-hic/GSE63525_HeLa_HiCCUPS_looplist.annotated.bed.gz,GSE63525_HeLa_HiCCUPS_looplist,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/insitu-hic/GSE63525_HMEC_HiCCUPS_looplist.annotated.bed.gz,GSE63525_HMEC_HiCCUPS_looplist,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/insitu-hic/GSE63525_HUVEC_HiCCUPS_looplist.annotated.bed.gz,GSE63525_HUVEC_HiCCUPS_looplist,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/insitu-hic/GSE63525_IMR90_HiCCUPS_looplist.annotated.bed.gz,GSE63525_IMR90_HiCCUPS_looplist,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/insitu-hic/GSE63525_K562_HiCCUPS_looplist.annotated.bed.gz,GSE63525_K562_HiCCUPS_looplist,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/insitu-hic/GSE63525_KBM7_HiCCUPS_looplist.annotated.bed.gz,GSE63525_KBM7_HiCCUPS_looplist,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/insitu-hic/GSE63525_NHEK_HiCCUPS_looplist.annotated.bed.gz,GSE63525_NHEK_HiCCUPS_looplist,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/insitu-hic/GSE63525_GM12878_100kb_inter_MAPQGE30_SQRTVC.bed.gz,GSE63525_GM12878_100kb_inter_MAPQGE30_SQRTVC,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/e-mtab-2323/TS5_CD34_promoter-other_significant_interactions.bed.gz,TS5_CD34_promoter_other_significant_interactions,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/e-mtab-2323/TS5_CD34_promoter-promoter_significant_interactions.bed.gz,TS5_CD34_promoter_promoter_significant_interactions,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/e-mtab-2323/TS5_GM12878_promoter-other_significant_interactions.bed.gz,TS5_GM12878_promoter_other_significant_interactions,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/e-mtab-2323/TS5_GM12878_promoter-promoter_significant_interactions.bed.gz,TS5_GM12878_promoter_promoter_significant_interactions,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/ucsc/encode/wgEncodeAwgDnaseMasterSites.bed.gz,wgEncodeAwgDnaseMasterSites,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/ucsc/encode/wgEncodeRegDnaseClusteredV3.bed.gz,wgEncodeRegDnaseClusteredV3,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/fantom5/fantom5PermissiveEnhancers.bed.gz,fantom5PermissiveEnhancers,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/dbsuper/all/allDbSuperEnhancers.bed.gz,allDbSuperEnhancers,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/dnase/roadmapDnaseDyadic.bed.gz,roadmapDnaseDyadic,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/dnase/roadmapDnaseEnh.bed.gz,roadmapDnaseEnh,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/dnase/roadmapDnaseProm.bed.gz,roadmapDnaseProm,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/ucsc/encode/wgEncodeRegTfbsClusteredWithCellsV3.bed.gz,wgEncodeRegTfbsClusteredWithCellsV3,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/ucsc/encode/wgEncodeAwgSegmentationCombinedGm12878.bed.gz,wgEncodeAwgSegmentationCombinedGm12878,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/ucsc/encode/wgEncodeAwgSegmentationCombinedH1hesc.bed.gz,wgEncodeAwgSegmentationCombinedH1hesc,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/ucsc/encode/wgEncodeAwgSegmentationCombinedHelas3.bed.gz,wgEncodeAwgSegmentationCombinedHelas3,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/ucsc/encode/wgEncodeAwgSegmentationCombinedHepg2.bed.gz,wgEncodeAwgSegmentationCombinedHepg2,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/ucsc/encode/wgEncodeAwgSegmentationCombinedHuvec.bed.gz,wgEncodeAwgSegmentationCombinedHuvec,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/ucsc/encode/wgEncodeAwgSegmentationCombinedK562.bed.gz,wgEncodeAwgSegmentationCombinedK562,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E001_15_coreMarks_mnemonics.bed.gz,E001_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E002_15_coreMarks_mnemonics.bed.gz,E002_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E003_15_coreMarks_mnemonics.bed.gz,E003_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E004_15_coreMarks_mnemonics.bed.gz,E004_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E005_15_coreMarks_mnemonics.bed.gz,E005_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E006_15_coreMarks_mnemonics.bed.gz,E006_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E007_15_coreMarks_mnemonics.bed.gz,E007_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E008_15_coreMarks_mnemonics.bed.gz,E008_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E009_15_coreMarks_mnemonics.bed.gz,E009_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E010_15_coreMarks_mnemonics.bed.gz,E010_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E011_15_coreMarks_mnemonics.bed.gz,E011_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E012_15_coreMarks_mnemonics.bed.gz,E012_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E013_15_coreMarks_mnemonics.bed.gz,E013_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E014_15_coreMarks_mnemonics.bed.gz,E014_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E015_15_coreMarks_mnemonics.bed.gz,E015_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E016_15_coreMarks_mnemonics.bed.gz,E016_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E017_15_coreMarks_mnemonics.bed.gz,E017_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E018_15_coreMarks_mnemonics.bed.gz,E018_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E019_15_coreMarks_mnemonics.bed.gz,E019_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E020_15_coreMarks_mnemonics.bed.gz,E020_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E021_15_coreMarks_mnemonics.bed.gz,E021_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E022_15_coreMarks_mnemonics.bed.gz,E022_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E023_15_coreMarks_mnemonics.bed.gz,E023_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E024_15_coreMarks_mnemonics.bed.gz,E024_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E025_15_coreMarks_mnemonics.bed.gz,E025_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E026_15_coreMarks_mnemonics.bed.gz,E026_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E027_15_coreMarks_mnemonics.bed.gz,E027_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E028_15_coreMarks_mnemonics.bed.gz,E028_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E029_15_coreMarks_mnemonics.bed.gz,E029_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E030_15_coreMarks_mnemonics.bed.gz,E030_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E031_15_coreMarks_mnemonics.bed.gz,E031_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E032_15_coreMarks_mnemonics.bed.gz,E032_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E033_15_coreMarks_mnemonics.bed.gz,E033_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E034_15_coreMarks_mnemonics.bed.gz,E034_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E035_15_coreMarks_mnemonics.bed.gz,E035_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E036_15_coreMarks_mnemonics.bed.gz,E036_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E037_15_coreMarks_mnemonics.bed.gz,E037_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E038_15_coreMarks_mnemonics.bed.gz,E038_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E039_15_coreMarks_mnemonics.bed.gz,E039_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E040_15_coreMarks_mnemonics.bed.gz,E040_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E041_15_coreMarks_mnemonics.bed.gz,E041_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E042_15_coreMarks_mnemonics.bed.gz,E042_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E043_15_coreMarks_mnemonics.bed.gz,E043_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E044_15_coreMarks_mnemonics.bed.gz,E044_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E045_15_coreMarks_mnemonics.bed.gz,E045_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E046_15_coreMarks_mnemonics.bed.gz,E046_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E047_15_coreMarks_mnemonics.bed.gz,E047_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E048_15_coreMarks_mnemonics.bed.gz,E048_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E049_15_coreMarks_mnemonics.bed.gz,E049_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E050_15_coreMarks_mnemonics.bed.gz,E050_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E051_15_coreMarks_mnemonics.bed.gz,E051_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E052_15_coreMarks_mnemonics.bed.gz,E052_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E053_15_coreMarks_mnemonics.bed.gz,E053_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E054_15_coreMarks_mnemonics.bed.gz,E054_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E055_15_coreMarks_mnemonics.bed.gz,E055_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E056_15_coreMarks_mnemonics.bed.gz,E056_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E057_15_coreMarks_mnemonics.bed.gz,E057_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E058_15_coreMarks_mnemonics.bed.gz,E058_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E059_15_coreMarks_mnemonics.bed.gz,E059_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E061_15_coreMarks_mnemonics.bed.gz,E061_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E062_15_coreMarks_mnemonics.bed.gz,E062_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E063_15_coreMarks_mnemonics.bed.gz,E063_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E065_15_coreMarks_mnemonics.bed.gz,E065_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E066_15_coreMarks_mnemonics.bed.gz,E066_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E067_15_coreMarks_mnemonics.bed.gz,E067_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E068_15_coreMarks_mnemonics.bed.gz,E068_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E069_15_coreMarks_mnemonics.bed.gz,E069_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E070_15_coreMarks_mnemonics.bed.gz,E070_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E071_15_coreMarks_mnemonics.bed.gz,E071_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E072_15_coreMarks_mnemonics.bed.gz,E072_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E073_15_coreMarks_mnemonics.bed.gz,E073_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E074_15_coreMarks_mnemonics.bed.gz,E074_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E075_15_coreMarks_mnemonics.bed.gz,E075_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E076_15_coreMarks_mnemonics.bed.gz,E076_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E077_15_coreMarks_mnemonics.bed.gz,E077_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E078_15_coreMarks_mnemonics.bed.gz,E078_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E079_15_coreMarks_mnemonics.bed.gz,E079_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E080_15_coreMarks_mnemonics.bed.gz,E080_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E081_15_coreMarks_mnemonics.bed.gz,E081_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E082_15_coreMarks_mnemonics.bed.gz,E082_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E083_15_coreMarks_mnemonics.bed.gz,E083_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E084_15_coreMarks_mnemonics.bed.gz,E084_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E085_15_coreMarks_mnemonics.bed.gz,E085_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E086_15_coreMarks_mnemonics.bed.gz,E086_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E087_15_coreMarks_mnemonics.bed.gz,E087_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E088_15_coreMarks_mnemonics.bed.gz,E088_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E089_15_coreMarks_mnemonics.bed.gz,E089_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E090_15_coreMarks_mnemonics.bed.gz,E090_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E091_15_coreMarks_mnemonics.bed.gz,E091_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E092_15_coreMarks_mnemonics.bed.gz,E092_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E093_15_coreMarks_mnemonics.bed.gz,E093_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E094_15_coreMarks_mnemonics.bed.gz,E094_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E095_15_coreMarks_mnemonics.bed.gz,E095_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E096_15_coreMarks_mnemonics.bed.gz,E096_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E097_15_coreMarks_mnemonics.bed.gz,E097_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E098_15_coreMarks_mnemonics.bed.gz,E098_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E099_15_coreMarks_mnemonics.bed.gz,E099_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E100_15_coreMarks_mnemonics.bed.gz,E100_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E101_15_coreMarks_mnemonics.bed.gz,E101_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E102_15_coreMarks_mnemonics.bed.gz,E102_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E103_15_coreMarks_mnemonics.bed.gz,E103_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E104_15_coreMarks_mnemonics.bed.gz,E104_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E105_15_coreMarks_mnemonics.bed.gz,E105_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E106_15_coreMarks_mnemonics.bed.gz,E106_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E107_15_coreMarks_mnemonics.bed.gz,E107_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E108_15_coreMarks_mnemonics.bed.gz,E108_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E109_15_coreMarks_mnemonics.bed.gz,E109_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E110_15_coreMarks_mnemonics.bed.gz,E110_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E111_15_coreMarks_mnemonics.bed.gz,E111_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E112_15_coreMarks_mnemonics.bed.gz,E112_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E113_15_coreMarks_mnemonics.bed.gz,E113_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E114_15_coreMarks_mnemonics.bed.gz,E114_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E115_15_coreMarks_mnemonics.bed.gz,E115_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E116_15_coreMarks_mnemonics.bed.gz,E116_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E117_15_coreMarks_mnemonics.bed.gz,E117_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E118_15_coreMarks_mnemonics.bed.gz,E118_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E119_15_coreMarks_mnemonics.bed.gz,E119_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E120_15_coreMarks_mnemonics.bed.gz,E120_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E121_15_coreMarks_mnemonics.bed.gz,E121_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E122_15_coreMarks_mnemonics.bed.gz,E122_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E123_15_coreMarks_mnemonics.bed.gz,E123_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E124_15_coreMarks_mnemonics.bed.gz,E124_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E125_15_coreMarks_mnemonics.bed.gz,E125_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E126_15_coreMarks_mnemonics.bed.gz,E126_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E127_15_coreMarks_mnemonics.bed.gz,E127_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E128_15_coreMarks_mnemonics.bed.gz,E128_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/roadmap/chromhmm/E129_15_coreMarks_mnemonics.bed.gz,E129_15_coreMarks,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/ucsc/database/all_est.bed.gz,all_est,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/ucsc/database/gwasCatalog.bed.gz,gwasCatalog,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/ucsc/database/targetScanS.bed.gz,targetScanS,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/ucsc/database/tfbsConsSites.bed.gz,tfbsConsSites,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/ucsc/database/wgRna.bed.gz,wgRna,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/ucsc/database/phastConsElements46way.bed.gz,phastConsElements46way,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/ucsc/database/phastConsElements100way.bed.gz,phastConsElements100way,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/ucsc/database/gap.bed.gz,gap,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/ucsc/database/genomicSuperDups.bed.gz,genomicSuperDups,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/ucsc/database/rmsk.bed.gz,rmsk,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/ucsc/database/nestedRepeats.bed.gz,nestedRepeats,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/ucsc/database/simpleRepeat.bed.gz,simpleRepeat,bed,overlap,0 \\
    --custom /n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot/ucsc/database/microsat.bed.gz,microsat,bed,overlap,0 \\
    --input_file #{hotspotChrVepSigFile} \\
    --output_file #{hotspotChrVepSigOutFile}
    CONTENT
    file.puts content
  end
  cmd = <<-CMD
bsub -g /hot/vep \\
    -q i2b2_12h -W 12:0 \\
    -n 10 -R "span[hosts=1]" \\
    -o #{hotspotChrVepSigLsfOutFile} \\
    sh #{hotspotChrVepSigLsfShFile}
  CMD
  system cmd
end
