#!/usr/bin/env ruby

require 'pathname'

numCores = 1
baseDir = Pathname.new "/n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver"
rScript = baseDir + "script/2.noncodiver-hotspot_association.R"
annotFiles = Pathname.glob(baseDir + "hotspot/chrs/*/*.annotated.[0-9][0-9][0-9][0-9][0-9].txt").sort
annotFiles.each do |annotFile|
  expFile = annotFile.sub_ext(".asstested.txt")
  lsfout = expFile.sub_ext(".txt.lsfout")
  next if lsfout.exist?
        #-q mcore -n #{numCores} -W 700:0\\
  cmd =<<-CMD
      bsub \\
        -g /nd/hotspot \\
        -q short -n #{numCores} -W 12:0 \\
        -R "span[hosts=1]" \\
        -o #{lsfout} \\
        xvfb-run -a /opt/R-3.1.2/bin/Rscript #{rScript} #{annotFile} #{numCores}
  CMD
  system cmd
end
