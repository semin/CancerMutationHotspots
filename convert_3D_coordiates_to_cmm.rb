#!/usr/bin/env ruby

#<marker_set name="marker set 1">
#<marker id="1" x="-6.1267" y="17.44" z="-3.1338"  radius="0.35217"/>
#<marker id="2" x="1.5395" y="16.277" z="-3.0339" r="0" g="1" b="1"
               #radius="0.5" note="An example note"/>
               #<link id1="2" id2="1" r="1" g="1" b="0" radius="0.17609"/>
               #</marker_set>

markerCnt = 0
puts "<marker_set name=\"marker set 1\">"
ARGF.readlines.each_with_index do |line, i|
  if !line.empty?
    x, y, z = line.chomp.split("\t")
    puts "<marker id=\"#{i+1}\" x=\"#{x}\" y=\"#{y}\" z=\"#{z}\" radius=\"0.5\"/>"
    markerCnt += 1
  end
end

(1..(markerCnt-1)).each do |i|
  puts "<link id1=\"#{i+1}\" id2=\"#{i}\" radius=\"0.3\"/>"
end

puts "</marker_set>"


