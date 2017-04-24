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
if ( -d jobs ) then
   rm -fR jobs
endif
mkdir jobs
mv *_x???.txt jobs
cd jobs
foreach file ( $files )
   set sampleName = `basename $file .txt`
   mkdir $sampleName
   mv $file $sampleName
   if ( -e ../json.txt ) then
      cp -p ../json.txt $sampleName
      cp -p ../collisions.csv $sampleName
      cp -p ../lumis.csv $sampleName
      cp -p ../lumispath.csv $sampleName
   endif
   cd $sampleName
   set exeName = "qcdRates_"$sampleName
   echo MssmHbbHltPaths $file > $exeName
   chmod u+x $exeName
   $HOME/bin/qsub.sh $exeName $sampleName
   sleep 5
   cd -
end
exit
