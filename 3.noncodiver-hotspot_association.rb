#!/usr/bin/env ruby

require 'pathname'
require 'parallel'

baseDir = Pathname.new "/n/data1/hms/dbmi/park/semin/BiO/Research/Hotspot"
rScript = baseDir + "script/3.noncodiver-hotspot_association.R"
annotFiles = Pathname.glob(baseDir + "hotspot/chrs/*/*.annotated.[0-9][0-9][0-9][0-9][0-9].txt").sort
Parallel.each(annotFiles, :in_threads => 10) do |annotFile|
  expFile = annotFile.sub_ext(".ass.txt")
  lsfout = expFile.sub_ext(".txt.lsfout")
  next if lsfout.exist?
  cmd =<<-CMD
      bsub \\
        -g /hot/ass \\
        -q i2b2_12h -W 12:0 \\
        -o #{lsfout} \\
        xvfb-run -a /opt/R-3.1.2/bin/Rscript #{rScript} #{annotFile}
  CMD
  system cmd
end
