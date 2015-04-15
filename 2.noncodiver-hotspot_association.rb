#!/usr/bin/env ruby

require 'pathname'

numCores = 10
baseDir = Pathname.new "/n/data1/hms/dbmi/park/semin/BiO/Research/NoncoDiver"
rScript = baseDir + "script/3.noncodiver-hotspot_impact_on_expression.R"
annotFiles = Pathname.glob(baseDir + "hotspot/*.annotated.txt").sort
annotFiles.each do |annotFile|
  expFile = annotFile.sub_ext(".exptested.txt")
  lsfout = expFile.sub_ext(".txt.lsfout")
  next if lsfout.exist?
  cmd =<<-CMD
      bsub \\
        -g /nd/hotspot \\
        -q i2b2_1d -n #{numCores} \\
        -R "span[hosts=1]" \\
        -o #{lsfout} \\
        xvfb-run -a /opt/R-3.1.2/bin/Rscript #{rScript} #{annotFile}
  CMD
  system cmd
end
