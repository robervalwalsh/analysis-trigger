#!/bin/csh -f

set sampleDir = $1
set nSplit = 10

if ( $#argv < 1 ) then
   echo Need to give sample name
   exit
endif

if ( $#argv == 2 ) then
   set nSplit = $2
endif

cd $sampleDir
rm -f *_x???.txt
set baseDir = `basename $sampleDir`
set txtFile = `/bin/ls *$baseDir*.txt`
split.csh $nSplit $txtFile
set files = `/bin/ls *_x???.txt` 
foreach file ( $files )
   set sampleName = `basename $file .txt`
   echo MssmHbbHltPaths $file > qcdRates
   chmod u+x qcdRates
   $HOME/bin/qsub.sh qcdRates $sampleName
   sleep 5
end
exit
