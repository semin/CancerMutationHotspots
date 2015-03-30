#!/usr/bin/env ruby

def sanitize_vcfs
  headers = []
  ARGF.each_line do |line|
    if line.start_with?("##")
      headers << line.chomp
    elsif line.start_with?("#CHROM")
      headers << '##INFO=<ID=SC,Number=A,Type=Integer,Description="Sample count in genotypes, for each ALT allele, in the same order as listed">'
      headers << '##INFO=<ID=SN,Number=1,Type=Integer,Description="Total number of samples with somatic mutations">'
      headers << '##INFO=<ID=NAC,Number=A,Type=Integer,Description="Updated allele count in genotypes, for each ALT allele, in the same order as listed">'
      puts headers.join("\n")
      print line
    else
      chrom, pos, id, ref, alt, qual, filter, info, format, *samples = line.chomp.split("\t")
      genotypes = samples.map { |s| s.split(":")[0] }
      #pn = samples.size
      scs = alt.split(",").each_with_index.map { |a, i| genotypes.select { |g| g.split("/").include?((i+1).to_s) }.size }
      sn = genotypes.select { |g| !g.start_with?(".") }.size
      nacs = alt.split(",").each_with_index.map { |a, i| genotypes.map { |g| g.split("/") }.flatten.count((i+1).to_s) }
      #info = [info, "PN=#{pn}", "SC=#{scs.join(',')}", "NAC=#{nacs.join(',')}"].join(";")
      ninfo = [info, "SC=#{scs.join(',')}", "SN=#{sn}", "NAC=#{nacs.join(',')}"].join(";")
      puts [chrom, pos, id, ref, alt, qual, filter, ninfo, format, samples].flatten.join("\t")
    end
  end
end

sanitize_vcfs

