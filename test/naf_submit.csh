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
set newdir = /nfs/dust/cms/user/walsh/naf/jobs/jobs.$$
#mkdir jobs
mkdir $newdir
ln -s $newdir jobs
mv *_x???.txt jobs
set thisdir = `pwd`
echo $thisdir
cd jobs
foreach file ( $files )
   set sampleName = `basename $file .txt`
   mkdir $sampleName
   mv $file $sampleName
   if ( -e $thisdir/json.txt ) then
      cp -p $thisdir/json.txt $sampleName
      cp -p $thisdir/collisions.csv $sampleName
      cp -p $thisdir/lumis.csv $sampleName
      cp -p $thisdir/lumispath.csv $sampleName
   endif
   cd $sampleName
   set exeName = "rates_"$sampleName
   echo MssmHbbHltPaths $file > $exeName
   chmod u+x $exeName
   $HOME/bin/qsub.sh $exeName $sampleName
#   sleep 5
   cd -
end
exit
