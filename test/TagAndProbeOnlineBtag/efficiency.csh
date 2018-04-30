#!/bin/csh -f

if ( $#argv < 1 ) then
   echo "Need to give directory containing the histograms"
   exit
endif

set mydir = $1

set eras = ( "2017C-v1"  "2017C-v2"  "2017C-v3"  "2017C"  "2017D"  "2017E-v1"  "2017E-v2"  "2017F"  "2017E"  "2017CD"  "2017CDE"  "2017CDEF"  "2017CDE-v1"  "2017E-v2F" )

foreach era ( $eras )

   root -l -b -q  Efficiency.C\(\"$mydir\",\"$era\"\)

end

root -l -b -q  EffPlots.C\(\"$mydir\"\)
