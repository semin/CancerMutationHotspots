#!/usr/bin/env ruby

require 'logger'
require 'pathname'
require 'parallel'

$logger           = Logger.new(STDOUT)
$logger.level     = Logger::DEBUG

$java7_bin        = Pathname.new("/opt/java/jdk1.7.0_71/bin/java")
$home_dir         = Pathname.new("/home/sl279")
$base_dir         = $home_dir + "BiO/Research/NoncoDiver"
$script_dir       = $base_dir + "script"
$vcf_dir          = $base_dir + "vcf"
$install_dir      = $home_dir + "BiO/Install"
$gatk3_bin        = $install_dir + "GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar"
$gatk_bundle_dir  = $install_dir + "GATK-bundle"
$refseq           = $gatk_bundle_dir + "human_g1k_v37_decoy.fasta"
$dbsnp_dir        = $base_dir + "dbsnp"
$hotspot_dir      = $base_dir + "hotspot"
$chrom_len_b37    = $gatk_bundle_dir + "hg19.genome"
$chrs             = (1..22).to_a << "X" << "Y"
$chrchrs          = $chrs.map { |c| "chr#{c}" }
$tmpDir           = Pathname.new("/home/sl279/BiO/Temp"); $tmpDir.mkpath

def submit(cmd)
  out = `#{cmd}`.chomp
  if out =~ /Job\s{1}<(\d+)>/
    puts "Job #{$1} submitted."
    return Integer($1)
  else
    abort out
  end
end

def annotate_vcfs_with_dbsnp_rsids
  dbsnp = $base_dir + "dbsnp/All.rsid.exp.vcf.gz"
  vcfs = Pathname.glob($vcf_dir + "*.vcf.gz").sort
  Parallel.each_with_index(vcfs, :in_threads => 4) do |vcf, vi|
    dbsnp_vcf = vcf.dirname + vcf.basename(".gz").sub_ext(".rsid.vcf.gz")
    lsfout = dbsnp_vcf.sub_ext(".gz.lsfout")
    cmd = <<-CMD
    bsub \\
      -g /nd/rsid \\
      -q short -W 12:0 \\
      -o #{lsfout} \\
      "zcat #{vcf} | vcf-annotate -a #{dbsnp} -d key=INFO,ID=DBSNP,Number=A,Type=String,Description=\\"dbSNP142 RSID\\" -c CHROM,POS,INFO/DBSNP,REF,ALT | bgzip -c > #{dbsnp_vcf} && tabix -p vcf #{dbsnp_vcf}"
    CMD
    submit cmd
  end
end

def annotate_vcfs_with_1000genomes_afs
  tgaf = $base_dir + "1000genomes/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.afs.tsv.exp.gz"
  vcfs = Pathname.glob($vcf_dir + "*.rsid.vcf.gz").sort
  Parallel.each_with_index(vcfs, :in_threads => 4) do |vcf, vi|
    tgALL_vcf = vcf.dirname + vcf.basename(".gz").sub_ext(".tgALL.vcf.gz")
    tgEAS_vcf = tgALL_vcf.dirname + tgALL_vcf.basename(".gz").sub_ext(".tgEAS.vcf.gz")
    tgEUR_vcf = tgALL_vcf.dirname + tgEAS_vcf.basename(".gz").sub_ext(".tgEUR.vcf.gz")
    tgAFR_vcf = tgALL_vcf.dirname + tgEUR_vcf.basename(".gz").sub_ext(".tgAFR.vcf.gz")
    tgAMR_vcf = tgALL_vcf.dirname + tgAFR_vcf.basename(".gz").sub_ext(".tgAMR.vcf.gz")
    tgSAS_vcf = tgALL_vcf.dirname + tgAMR_vcf.basename(".gz").sub_ext(".tgSAS.vcf.gz")
    final_vcf = vcf.dirname + vcf.basename(".gz").sub_ext(".tg.vcf.gz")
    lsfout = final_vcf.sub_ext(".gz.lsfout")
    cmd = <<-CMD
    bsub \\
      -g /nd/tg \\
      -q short -W 12:0 \\
      -o #{lsfout} "
      zcat #{vcf} | vcf-annotate -a #{tgaf} -d key=INFO,ID=TG_ALL_AF,Number=A,Type=Float,Description=\\"Global Allele Frequency from 1000 Genomes Project\\" -c CHROM,POS,REF,ALT,INFO/TG_ALL_AF,-,-,-,-,- | bgzip -c > #{tgALL_vcf} &&
      tabix -p vcf #{tgALL_vcf} &&
      zcat #{tgALL_vcf} | vcf-annotate -a #{tgaf} -d key=INFO,ID=TG_EAS_AF,Number=A,Type=Float,Description=\\"EAS Allele Frequency from 1000 Genomes Project\\" -c CHROM,POS,REF,ALT,-,INFO/TG_EAS_AF,-,-,-,- | bgzip -c > #{tgEAS_vcf} &&
      tabix -p vcf #{tgEAS_vcf} &&
      zcat #{tgEAS_vcf} | vcf-annotate -a #{tgaf} -d key=INFO,ID=TG_EUR_AF,Number=A,Type=Float,Description=\\"EUR Allele Frequency from 1000 Genomes Project\\" -c CHROM,POS,REF,ALT,-,-,INFO/TG_EUR_AF,-,-,- | bgzip -c > #{tgEUR_vcf} &&
      tabix -p vcf #{tgEUR_vcf} &&
      zcat #{tgEUR_vcf} | vcf-annotate -a #{tgaf} -d key=INFO,ID=TG_AFR_AF,Number=A,Type=Float,Description=\\"AFR Allele Frequency from 1000 Genomes Project\\" -c CHROM,POS,REF,ALT,-,-,-,INFO/TG_AFR_AF,-,- | bgzip -c > #{tgAFR_vcf} &&
      tabix -p vcf #{tgAFR_vcf} &&
      zcat #{tgAFR_vcf} | vcf-annotate -a #{tgaf} -d key=INFO,ID=TG_AMR_AF,Number=A,Type=Float,Description=\\"AMR Allele Frequency from 1000 Genomes Project\\" -c CHROM,POS,REF,ALT,-,-,-,-,INFO/TG_AMR_AF,- | bgzip -c > #{tgAMR_vcf} &&
      tabix -p vcf #{tgAMR_vcf} &&
      zcat #{tgAMR_vcf} | vcf-annotate -a #{tgaf} -d key=INFO,ID=TG_SAS_AF,Number=A,Type=Float,Description=\\"SAS Allele Frequency from 1000 Genomes Project\\" -c CHROM,POS,REF,ALT,-,-,-,-,-,INFO/TG_SAS_AF | bgzip -c > #{tgSAS_vcf} &&
      rm -rf #{tgALL_vcf} #{tgEAS_vcf} #{tgEUR_vcf} #{tgAFR_vcf} #{tgAMR_vcf} &&
      rm -rf #{tgALL_vcf}.tbi #{tgEAS_vcf}.tbi #{tgEUR_vcf}.tbi #{tgAFR_vcf}.tbi #{tgAMR_vcf}.tbi &&
      mv #{tgSAS_vcf} #{final_vcf} &&
      tabix -p vcf #{final_vcf}
      "
    CMD
    submit cmd
  end
end

def annotate_vcfs_with_esp6500si_afs
  espAfFile = $base_dir + "nhlbi" + "ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.afs.sorted.txt.gz"
  vcfs = Pathname.glob($vcf_dir + "*.tg.vcf.gz").sort
  Parallel.each_with_index(vcfs, :in_threads => 4) do |vcf, vi|
    espTA_vcf = vcf.dirname + vcf.basename(".gz").sub_ext(".espTA.vcf.gz")
    espEA_vcf = espTA_vcf.dirname + espTA_vcf.basename(".gz").sub_ext(".espEA.vcf.gz")
    espAA_vcf = espTA_vcf.dirname + espEA_vcf.basename(".gz").sub_ext(".espAA.vcf.gz")
    final_vcf = vcf.dirname + vcf.basename(".gz").sub_ext(".esp.vcf.gz")
    lsfout = final_vcf.sub_ext(".gz.lsfout")
    cmd = <<-CMD
    bsub \\
      -g /nd/esp \\
      -q short -W 12:0 \\
      -o #{lsfout} "
      zcat #{vcf} | vcf-annotate -a #{espAfFile} -d key=INFO,ID=TA_AF,Number=A,Type=Float,Description=\\"Total American Allele Frequency from NHLBI ESP6500SI-V2-SSA137\\" -c CHROM,POS,REF,ALT,-,-,INFO/TA_AF | bgzip -c > #{espTA_vcf} &&
      tabix -p vcf #{espTA_vcf} &&
      zcat #{espTA_vcf} | vcf-annotate -a #{espAfFile} -d key=INFO,ID=EA_AF,Number=A,Type=Float,Description=\\"European American Allele Frequency from NHLBI ESP6500SI-V2-SSA137\\" -c CHROM,POS,REF,ALT,INFO/EA_AF,-,- | bgzip -c > #{espEA_vcf} &&
      tabix -p vcf #{espEA_vcf} &&
      zcat #{espEA_vcf} | vcf-annotate -a #{espAfFile} -d key=INFO,ID=AA_AF,Number=A,Type=Float,Description=\\"African American Allele Frequency from NHLBI ESP6500SI-V2-SSA137\\" -c CHROM,POS,REF,ALT,-,INFO/AA_AF,- | bgzip -c > #{espAA_vcf} &&
      rm -rf #{espTA_vcf} #{espEA_vcf} &&
      rm -rf #{espTA_vcf}.tbi #{espEA_vcf}.tbi &&
      mv #{espAA_vcf} #{final_vcf} &&
      tabix -p vcf #{final_vcf}
      "
    CMD
    submit cmd
  end
end

def annotate_snp_vcfs_with_cadd
  cadd = $base_dir + "cadd/v1.2/whole_genome_SNVs.tsv.gz"
  vcfs = Pathname.glob($vcf_dir + "cancer/*/*.esp.vcf.gz").sort
  Parallel.each_with_index(vcfs, :in_threads => 4) do |vcf, vi|
    cadd_raw_vcf = vcf.dirname + vcf.basename(".gz").sub_ext(".caddRaw.vcf.gz")
    cadd_phred_vcf = vcf.dirname + cadd_raw_vcf.basename(".gz").sub_ext(".caddPhred.vcf.gz")
    final_vcf = vcf.dirname + vcf.basename(".gz").sub_ext(".cadd.vcf.gz")
    lsfout = final_vcf.sub_ext(".gz.lsfout")
      #-q i2b2_7d \\
      #-q long -W 700:0 \\
    cmd = <<-CMD
    bsub \\
      -g /nd/cadd \\
      -q park_7d \\
      -o #{lsfout} \\
      "zcat #{vcf} | vcf-annotate -a #{cadd} -d key=INFO,ID=CADD_RAW,Number=A,Type=Float,Description=\\"CADD Raw Score\\" -c CHROM,POS,REF,ALT,INFO/CADD_RAW,- | bgzip -c > #{cadd_raw_vcf} &&
      tabix -p vcf #{cadd_raw_vcf} &&
      zcat #{cadd_raw_vcf} | vcf-annotate -a #{cadd} -d key=INFO,ID=CADD_PHRED,Number=A,Type=Float,Description=\\"CADD Phred Score\\" -c CHROM,POS,REF,ALT,-,INFO/CADD_PHRED | bgzip -c > #{cadd_phred_vcf}
      rm -rf #{cadd_raw_vcf} #{cadd_raw_vcf}.tbi &&
      mv #{cadd_phred_vcf} #{final_vcf} &&
      tabix -p vcf #{final_vcf}"
    CMD
    submit cmd
  end
end

def annotate_vcfs_with_vep
  vep_bin = Pathname.new "/home/sl279/vep/variant_effect_predictor.pl"
  vep_dir = Pathname.new "/home/sl279/.vep"
  encode_dcc_dir = $base_dir + "ucsc/encode"
  ucsc_db_dir = $base_dir + "ucsc/database"
  fantom5_dir = $base_dir + "fantom5"
  dbsuper_dir = $base_dir + "dbsuper/all"
  vcfs = Pathname.glob($vcf_dir + "cancer/*/*.cadd.vcf.gz").sort
  Parallel.each_with_index(vcfs, :in_threads => 4) do |vcf, vi|
    vep_vcf = vcf.dirname + vcf.basename(".gz").sub_ext(".vep.vcf.gz")
    lsfout = vep_vcf.sub_ext(".gz.lsfout")
    #next if lsfout.exist?
    cmd = <<-CMD
    bsub \\
      -g /nd/cancer/vep \\
      -q i2b2_7d \\
      -n 8 -R "span[hosts=1]" \\
      -o #{lsfout} \\
      "perl #{vep_bin} \\
        --force_overwrite \\
        --fork 8 \\
        --offline \\
        --no_stats \\
        --everything \\
        --check_existing \\
        --total_length \\
        --allele_number \\
        --no_escape \\
        --gencode_basic \\
        --vcf \\
        --assembly GRCh37 \\
        --dir #{vep_dir} \\
        --fasta #{vep_dir}/homo_sapiens/76_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \\
        --custom #{ucsc_db_dir}/gap.bed.gz,gap,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/genomicSuperDups.bed.gz,genomicSuperDups,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/gwasCatalog.bed.gz,gwasCatalog,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/targetScanS.bed.gz,targetScanS,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/tfbsConsSites.bed.gz,tfbsConsSites,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/wgRna.bed.gz,wgRna,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/phastConsElements46way.bed.gz,phastConsElements46way,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/phastConsElements100way.bed.gz,phastConsElements100way,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/rmsk.bed.gz,rmsk,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/nestedRepeats.bed.gz,nestedRepeats,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/simpleRepeat.bed.gz,simpleRepeat,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/microsat.bed.gz,microsat,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/wgEncodeAwgDnaseMasterSites.bed.gz,wgEncodeAwgDnaseMasterSites,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/wgEncodeRegDnaseClusteredV3.bed.gz,wgEncodeRegDnaseClusteredV3,bed,overlap,0 \\
        --custom #{fantom5_dir}/fantom5PermissiveEnhancers.bed.gz,fantom5PermissiveEnhancers,bed,overlap,0 \\
        --custom #{dbsuper_dir}/allDbSuperEnhancers.bed.gz,allDbSuperEnhancers,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/wgEncodeRegTfbsClusteredWithCellsV3.bed.gz,wgEncodeRegTfbsClusteredWithCellsV3,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/wgEncodeAwgSegmentationCombinedGm12878.bed.gz,wgEncodeAwgSegmentationCombinedGm12878,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/wgEncodeAwgSegmentationCombinedH1hesc.bed.gz,wgEncodeAwgSegmentationChromhmmH1hesc,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/wgEncodeAwgSegmentationCombinedHelas3.bed.gz,wgEncodeAwgSegmentationCombinedHelas3,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/wgEncodeAwgSegmentationCombinedHepg2.bed.gz,wgEncodeAwgSegmentationCombinedHepg2,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/wgEncodeAwgSegmentationCombinedHuvec.bed.gz,wgEncodeAwgSegmentationCombinedHuvec,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/wgEncodeAwgSegmentationCombinedK562.bed.gz,wgEncodeAwgSegmentationCombinedK562,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/distalDhsToPromoterDhs.bed.gz,distalDhsToPromoterDhs,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/dhsToGeneExpression.bed.gz,dhsToGeneExpression,bed,overlap,0 \\
        --custom #{fantom5_dir}/fantom5EnhancerTssAssociations.bed.gz,fantom5EnhancerTssAssociations,bed,overlap,0 \\
        --custom #{dbsuper_dir}/allDbSuperEnhancerGeneAssociations.bed.gz,allDbSuperEnhancerGeneAssociations,bed,overlap,0 \\
        --input_file #{vcf} \\
        --output_file STDOUT | bgzip -c > #{vep_vcf} && tabix -p vcf #{vep_vcf}"
    CMD
    submit cmd
  end
end

def split_vcfs_by_chr
  vcfs = Pathname.glob($vcf_dir + "*/*.vep.vcf.gz").sort
  vcfs.each do |vcf|
    $chrs.each do |chr|
      vcf_chr_dir = vcf.dirname + "chrs"; vcf_chr_dir.mkpath
      vcf_chr_file = vcf_chr_dir + vcf.basename(".gz").sub_ext(".chr#{chr}.vcf.gz")
      lsfout = vcf_chr_file.sub_ext(".gz.lsfout")
      cmd = <<-CMD
          bsub \\
            -g /nd/split \\
            -q short -W 12:0 \\
            -R "rusage[mem=10000]" -M 10000000 \\
            -o #{lsfout} \\
              #{$java7_bin} -Xmx10g -jar #{$gatk3_bin} \\
              -T SelectVariants \\
              -R #{$refseq} \\
              --variant #{vcf} \\
              -o #{vcf_chr_file} \\
              -L #{chr}
      CMD
      submit cmd
    end
  end
end

def sanitize_vcfs(vcfs)
  vcfs.each do |vcf|
    svcf = vcf.dirname + vcf.basename(".gz").sub_ext(".sanitized.vcf.gz")
    lsfout = svcf.sub_ext(".gz.lsfout")
    cmd = <<-CMD
    bsub \\
      -g /nd/sani \\
      -q short -W 12:0 \\
      -o #{lsfout} \\
      "zcat #{vcf} | ruby /n/data1/hms/dbmi/park/semin/BiO/Research/Recuronco/script/sanitize_vcfs.rb | bgzip -c > #{svcf} && tabix -p vcf #{svcf}"
    CMD
    submit cmd
  end
end

def split_vcfs_by_recurrence
  vcfs = Pathname.glob($vcf_dir + "pancan/*.sanitized.vcf.gz").sort
  vcfs.each do |vcf|
    puts vcf
    #rvcf = vcf.dirname + vcf.basename(".gz").sub_ext(".recurrent.vcf")
    ovcf = vcf.dirname + vcf.basename(".gz").sub_ext(".oneoff.vcf")
    #rvcf.open('w') do |rfile|
      ovcf.open('w') do |ofile|
        `zcat #{vcf}`.split("\n").each do |line|
          if line.start_with?("#")
            #rfile.puts line
            ofile.puts line
          else
            sn = line.match(/SN=(\d+)/)[1].to_i
            #rfile.puts line if sn > 1
            ofile.puts line if sn == 1
          end
        end
        #system "bgzip #{rvcf}; tabix -p vcf #{rvcf}.gz"
        system "bgzip #{ovcf}; tabix -p vcf #{ovcf}.gz"
      end
    end
  #end
end

def merge_chr_vcfs
  $chrs.each do |chr|
    vcfs = Pathname.glob($vcf_dir + "CANCER/*/chrs/*.vep.chr#{chr}.vcf.gz").sort
    ovcf = $vcf_dir + "pancan/TCGA_16_Cancer_Types.wgs.somatic.chr#{chr}.vcf.gz"
    lsfout = ovcf.sub_ext(".gz.lsfout")
    cmd = <<-CMD
    bsub \\
      -g /nd/merge \\
      -q short -W 12:0 \\
      -o #{lsfout} \\
      "vcf-merge #{vcfs.join(" ")} | bgzip -c > #{ovcf} && tabix -p vcf #{ovcf}"
    CMD
    submit cmd
  end
end

def extract_mutation_count
  vcfs = Pathname.glob($vcf_dir + "pancan/*sanitized.vcf.gz").sort
  Parallel.each(vcfs, :in_processes => 8) do |vcf|
    mcntFile = $hotspot_dir + "TCGA" + vcf.basename(".gz").sub_ext(".scnt.alt_expanded.txt")
    puts mcntFile
    mcntFile.open('w') do |file|
      sampleIds = []
      `zcat #{vcf}`.split("\n").each do |line|
        if line.start_with?("#CHROM")
          #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  H_LS-A1-A0SM-01A-11D-A19H-09
          sampleIds = line.chomp.split("\t")[9..-1]
        elsif line.start_with?("#")
          next
        else
          chrom, pos, id, ref, alt, qual, filter, info, format, *genos = line.chomp.split("\t")
          alts = alt.split(",")
          genos = genos.map { |g| g.split(":")[0].split("/") }
          caddrs = info.match(/CADD_RAW=(\S+?);/)[1].split(",")
          caddps = info.match(/CADD_PHRED=(\S+?);/)[1].split(",")
          mcaddrs = (alts.size != caddrs.size) ? Array.new(alts.size, caddrs[0]) : caddrs
          mcaddps = (alts.size != caddps.size) ? Array.new(alts.size, caddps[0]) : caddps
          nalts = []
          scnts = []
          nsids = []
          ncaddr = []
          ncaddp = []
          alts.each_with_index do |a, ai|
            gt = (ai + 1).to_s
            sids = sampleIds.values_at(*(genos.each_index.select{ |i| genos[i].include?(gt) }))
            if (sids.size > 0)
              file.puts [chrom, pos, ref, a, sids.join(",")].join("\t")
            end
            #if (sids.size > 0)
              #nalts << a
              #nsids[ai] = sids.join(";")
              #scnts[ai] = sids.size
              #ncaddr[ai] = mcaddrs[ai]
              #ncaddp[ai] = mcaddps[ai]
            #end
          end
          #file.puts [chrom, pos, ref,
                     #nalts.join(","), scnts.join(","), nsids.join(","),
                     #ncaddr.join(","), ncaddp.join(",")].join("\t")
        end
      end
    end
  end
end

def expand_vep_out
  vepOutFiles = Pathname.glob($vcf_dir + "pancan/*.sig2.vep_out.txt")
  #vepOutFiles.each do |vepOutFile|
  Parallel.each(vepOutFiles, :in_processes => 8) do |vepOutFile|
    puts vepOutFile
    vepExpOutFile = vepOutFile.sub_ext(".exp.txt")
    vepExpOutFile.open('w') do |file|
      extKeys = []
      vepOutFile.each_line do |line|
        if line.start_with?("##")
          if line =~ /^##\s(\S+)\s:/
            extKeys << $1
          end
        elsif line.start_with?("#Uploaded_variation")
          oriHeaders = line.chomp.split("\t")
          oriHeaders[0] = "Uploaded_variation"
          newHeaders = [oriHeaders[0..-2], extKeys].flatten
          file.puts newHeaders.join("\t")
        else
          extVals = Hash.new("")
          oriCols = line.chomp.split("\t")
          oriExts = oriCols[-1].split(";")
          oriExts.each do |oriExt|
            oriExtKey, oriExtVal = oriExt.split("=")
            extVals[oriExtKey] = oriExtVal
          end
          newCols = [oriCols[0..-2], extKeys.map { |k| extVals[k] }].flatten
          file.puts newCols.join("\t")
        end
      end
    end
  end
end

def extract_wgs_mutations_from_stratton_call_sets
  selectedSamples = []
  summaryFile = $base_dir + "sanger/samples_summary.txt"
  summaryFile.each_line do |line|
    cancerType, sampleName, seqType, dataSource = line.chomp.split("\t")
    if seqType == "Whole genome" && !dataSource.match(/TCGA/)
      selectedSamples << sampleName
    end
  end

  mutationCallFiles = Pathname.glob($base_dir + "sanger/somatic_mutation_data/*/*clean_somatic_mutations_for_signature_analysis.txt")
  mutationCallFiles.each do |mutationCallFile|
    mutationCallWgsFile = mutationCallFile.sub_ext(".wgs.txt")
    puts mutationCallWgsFile
    mutationCallWgsFile.open('w') do |file|
      mutationCallFile.each_line do |line|
        sampleName, mutType, chrom, start, stop, ref, alt = line.chomp.split("\t")
        file.print line if (selectedSamples.include?(sampleName) && mutType == "subs")
      end
    end
  end
end

if __FILE__ == $0
  #annotate_vcfs_with_dbsnp_rsids
  #annotate_vcfs_with_1000genomes_afs
  #annotate_vcfs_with_esp6500si_afs
  #annotate_snp_vcfs_with_cadd
  #annotate_vcfs_with_vep
  #split_vcfs_by_chr
  #sanitize_vcfs(Pathname.glob($vcf_dir + "pancan/*.somatic.chr*.vcf.gz").sort)
  #split_vcfs_by_recurrence
  #extract_mutation_count
  #expand_vep_out
  #extract_wgs_mutations_from_stratton_call_sets
end
