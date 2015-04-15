#!/usr/bin/env ruby

require 'pathname'
require 'parallel'

rscript = Pathname.new("/home/sl279/BiO/Research/NoncoDiver/script/postprocess_insitu-hic_data.R")
obsFiles = Pathname.glob("/home/sl279/BiO/Research/NoncoDiver/insitu-hic/*/100kb*/*/MAPQGE30/*.RAWobserved")
Parallel.each(obsFiles, :in_threads => 8) do |obsFile|
  #"chr1_10kb.RAWobserved"
  #"chr10_11_100kb.RAWobserved"
  chr, binPrefix, binPostfix = $1, $2, $3 if obsFile.basename.to_s =~ /(.*)\_(\d+)(\w+)b\.RAWobserved/
  chrs = chr.split("_")
  chrs = case chrs.size
         when 1
           [chrs[0], chrs[0]]
         when 2
           [chrs[0], "chr#{chrs[1]}"]
         else
           exit
         end
  chrA = chrs[0]
  chrB = chrs[1]
  binPostfixNum = case binPostfix
                  when "k"
                    1000
                  when "m"
                    1000000
                  end
  binSize = binPrefix.to_i * binPostfixNum
  normTypes = %w[KR SQRTVC VC]
  normTypes.each do |normType|
    normFileA = obsFile.dirname + "#{chrA}_#{binPrefix}#{binPostfix}b.#{normType}norm"
    normFileB = obsFile.dirname + "#{chrB}_#{binPrefix}#{binPostfix}b.#{normType}norm"
    if (normFileA.size == 0)
      warn "#{normFileA}'s size is 0!"
      next
    end
    if (normFileB.size == 0)
      warn "#{normFileB}'s size is 0!"
      next
    end
    expFile = obsFile.sub_ext("#{normType}expected")
    expFile = "NA" unless expFile.exist?
    outFile = obsFile.sub_ext(".#{normType}.bed")
    lsfout = outFile.sub_ext(".bed.lsfout")
    next if lsfout.exist?
    cmd = <<-CMD
bsub \\
  -g /gcc/germ/hic \\
  -q short -W 12:0 \\
  -o #{lsfout} \\
  /opt/R-3.1.2/bin/Rscript #{rscript} #{binSize} #{obsFile} #{normFileA} #{normFileB} #{expFile} #{outFile}
    CMD
    system cmd
  end
end
