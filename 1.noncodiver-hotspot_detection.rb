#!/usr/bin/env ruby

require 'pathname'

distCut = 50
hotspotMargin = 25
poisMargin = 5000
numCores = 10
baseDir = Pathname.new "/n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver"
rScript = baseDir + "script/1.noncodiver-hotspot_detection.R"
scntFiles = Pathname.glob(baseDir + "hotspot/*.scnt.txt").sort
scntFiles.each do |scntFile|
  hotspotFile = scntFile.sub_ext(".hotspot#{distCut}.txt")
  hotspotVepFile = hotspotFile.sub_ext(".vep_in.txt")
  lsfout = hotspotFile.sub_ext(".txt.lsfout")
  next if lsfout.exist?
  cmd =<<-CMD
      bsub \\
        -g /nd/hotspot \\
        -q i2b2_1d -n #{numCores} \\
        -R "span[hosts=1]" \\
        -o #{lsfout} \\
        xvfb-run -a /opt/R-3.1.2/bin/Rscript #{rScript} #{scntFile} #{hotspotFile} #{hotspotVepFile} #{distCut} #{hotspotMargin} #{poisMargin} #{numCores}
  CMD
  system cmd
end
