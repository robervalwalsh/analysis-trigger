#!/bin/csh -f

set pthat = $1
set nSplit = 10

if ( $#argv < 1 ) then
   echo Need to give sample name
   exit
endif

if ( $#argv == 2 ) then
   set nSplit = $2
endif

cd $pthat
rm -f *_x???.txt
split.csh $nSplit $pthat.txt
set files = `/bin/ls *_x???.txt` 
foreach file ( $files )
   set sampleName = `basename $file .txt`
   echo MssmHbbHltPaths $file > qcdRates
   chmod u+x qcdRates
   ../qsub.sh qcdRates $sampleName
   sleep 5
end
exit
