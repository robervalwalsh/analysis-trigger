#!/bin/csh -f

if ( $#argv < 5 ) then
   echo "Need to give era and number of splits, btag algo, tag and probe btag WPs, optionally the reco type"
   exit
endif

set era = $1
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
set mydir = $PWD/$reco/$name
if ( ! -d $mydir ) then
   mkdir -p $mydir
endif

cd $mydir

set json = "$mytestdir/certifiedlumis/$reco/Cert_"$era".json"
set ntuples = "$mytestdir/ntupleslists/$reco/BTagCSV_"$era".txt"

if ( ! -e $json ) then
   echo $json does not exist
   exit
endif

if ( ! -e $ntuples ) then
   echo $ntuples does not exist
   exit
endif

set ptmin =  (40  160 220 400)
set ptmax =  (160 220 400 1000)

set hlt    = ( "HLT_DoublePFJets40_CaloBTagCSV_p33_v" "HLT_DoublePFJets100_CaloBTagCSV_p33_v" "HLT_DoublePFJets200_CaloBTagCSV_p33_v" "HLT_DoublePFJets350_CaloBTagCSV_p33_v" )
set ol1    = ( "hltL1DoubleJet40er3p0" "hltL1DoubleJet100er3p0" "hltL1DoubleJet120er3p0" "hltL1DoubleJet120er3p0" )  
set ocalo  = ( "hltDoubleCaloBJets30eta2p3" "hltDoubleCaloBJets100eta2p3" "hltDoubleCaloBJets100eta2p3" "hltDoubleCaloBJets100eta2p3" )
set obtag  = ( "hltBTagCalo30x8CSVp0p92SingleWithMatching" "hltBTagCalo80x6CSVp0p92SingleWithMatching" "hltBTagCalo80x6CSVp0p92SingleWithMatching" "hltBTagCalo80x6CSVp0p92SingleWithMatching" )
set opf    = ( "hltDoublePFJets40Eta2p3" "hltDoublePFJets100Eta2p3" "hltDoublePFJets200Eta2p3" "hltDoublePFJets350Eta2p3" )

# change in L1 seed
if ( $era == "2017C-v1" || $era == "2017C-v2" ) then
   set ol1[3] = "hltL1DoubleJet112er3p0"
   set ol1[4] = "hltL1DoubleJet112er3p0"
endif


@ count = 1
while ( $count < $#ptmin + 1 )
   set cfgfile = "tep_"$era"_"$ptmin[$count]"to"$ptmax[$count]".cfg"
   set outfile = "histograms_"$era"_tep_"$ptmin[$count]"to"$ptmax[$count]".root"
  
### CONFIG FILE

   echo $cfgfile

   cp -pf $mypwd/tag_and_probe.cfg $cfgfile
   ln -s $json .
   ln -s $ntuples .
   set localjson = `basename $json`
   set localntpl = `basename $ntuples`
   
   sed -i -e "s/ERA/$era/g" $cfgfile
   sed -i -e "s/PTMIN/$ptmin[$count]/g" $cfgfile
   sed -i -e "s/PTMAX/$ptmax[$count]/g" $cfgfile
   sed -i -e "s/HLT/$hlt[$count]/g" $cfgfile
   sed -i -e "s/OL1/$ol1[$count]/g" $cfgfile
   sed -i -e "s/OCALO/$ocalo[$count]/g" $cfgfile
   sed -i -e "s/OBTAG/$obtag[$count]/g" $cfgfile
   sed -i -e "s/OPF/$opf[$count]/g" $cfgfile
   sed -i -e "s/JSON/$localjson/g" $cfgfile
   sed -i -e "s/TAGWP/$wptag/g" $cfgfile
   sed -i -e "s/PROBEWP/$wpprobe/g" $cfgfile
   sed -i -e "s/ALGBTAG/$btagalgo/g" $cfgfile
      
   naf_submit.csh $localntpl TagAndProbeOnlineBtag $cfgfile $splits $localjson
   
   rm -f $mydir/$cfgfile $mydir/$localjson $mydir/$localntpl

   @ count ++
end
