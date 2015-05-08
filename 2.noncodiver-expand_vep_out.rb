#!/usr/bin/env ruby

extKeys = []

ARGF.each_line do |line|
  if line.start_with?("##")
    if line =~ /^##\s(\S+)\s:/
      extKeys << $1
    end
  elsif line.start_with?("#Uploaded_variation")
    oriHeaders = line.chomp.split("\t")
    oriHeaders[0] = "Uploaded_variation"
    newHeaders = [oriHeaders[0..-2], extKeys].flatten
    puts newHeaders.join("\t")
  else
    extVals = Hash.new("")
    oriCols = line.chomp.split("\t")
    oriExts = oriCols[-1].split(";")
    oriExts.each do |oriExt|
      oriExtKey, oriExtVal = oriExt.split("=")
      extVals[oriExtKey] = oriExtVal
    end
    newCols = [oriCols[0..-2], extKeys.map { |k| extVals[k] }].flatten
    puts newCols.join("\t")
  end
end
