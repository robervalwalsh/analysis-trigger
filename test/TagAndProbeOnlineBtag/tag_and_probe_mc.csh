#!/bin/csh -f

if ( $#argv < 5 ) then
   echo "Need to give process and number of splits, btag algo, tag and probe btag WPs, optionally the reco type"
   exit
endif

set process = $1
set splits = $2
set btagalgo = $3
set wptag = $4
set wpprobe = $5

set name = $btagalgo"_"$wptag"_"$wpprobe

set reco = "ReReco17Nov2017_2017"
if ( $#argv == 6 ) then
   set reco = $6
endif

set mypwd = $PWD
set mytestdir = $PWD/..
set mydir = $PWD/MonteCarlo/$name
if ( ! -d $mydir ) then
   mkdir -p $mydir
endif

cd $mydir

set ntuples = "$mytestdir/ntupleslists/MonteCarlo/$process.txt"

if ( ! -e $ntuples ) then
   echo $ntuples does not exist
   exit
endif

set ptmin =  (40)
set ptmax =  (1000)

set hlt    = ( "HLT_DoublePFJets40_CaloBTagCSV_p33_v" )
set ol1    = ( "hltL1DoubleJet40er3p0" )  
set ocalo  = ( "hltDoubleCaloBJets30eta2p3" )
set obtag  = ( "hltBTagCalo30x8CSVp0p92SingleWithMatching" )
set opf    = ( "hltDoublePFJets40Eta2p3" )

@ count = 1
while ( $count < $#ptmin + 1 )
   set cfgfile = "tep_"$process"_"$ptmin[$count]"to"$ptmax[$count]".cfg"
   set outfile = "histograms_"$process"_tep_"$ptmin[$count]"to"$ptmax[$count]".root"
  
### CONFIG FILE

   echo $cfgfile

   cp -pf $mypwd/tag_and_probe_mc.cfg $cfgfile
   ln -s $ntuples .
   set localntpl = `basename $ntuples`
   
   sed -i -e "s/PROCESS/$process/g" $cfgfile
   sed -i -e "s/PTMIN/$ptmin[$count]/g" $cfgfile
   sed -i -e "s/PTMAX/$ptmax[$count]/g" $cfgfile
   sed -i -e "s/HLT/$hlt[$count]/g" $cfgfile
   sed -i -e "s/OL1/$ol1[$count]/g" $cfgfile
   sed -i -e "s/OCALO/$ocalo[$count]/g" $cfgfile
   sed -i -e "s/OBTAG/$obtag[$count]/g" $cfgfile
   sed -i -e "s/OPF/$opf[$count]/g" $cfgfile
   sed -i -e "s/TAGWP/$wptag/g" $cfgfile
   sed -i -e "s/PROBEWP/$wpprobe/g" $cfgfile
   sed -i -e "s/ALGBTAG/$btagalgo/g" $cfgfile
      
   naf_submit.csh $localntpl TagAndProbeOnlineBtag $cfgfile $splits
   
   rm -f $mydir/$cfgfile $mydir/$localntpl

   @ count ++
end
