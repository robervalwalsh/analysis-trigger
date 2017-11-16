#!/bin/csh -f

if ( $#argv < 3 ) then
   echo Need to give era, json and number of splits
   exit
endif

set era = $1
set json = $2
set splits = $3

set ptmin =  (40  160 220 400)
set ptmax =  (160 220 400 1000)

set hlt    = ( "HLT_DoublePFJets40_CaloBTagCSV_p33_v" "HLT_DoublePFJets100_CaloBTagCSV_p33_v" "HLT_DoublePFJets200_CaloBTagCSV_p33_v" "HLT_DoublePFJets350_CaloBTagCSV_p33_v" )
# for runs 2017 C-v3,D,E,F (below)
set ol1    = ( "hltL1DoubleJet40er3p0" "hltL1DoubleJet100er3p0" "hltL1DoubleJet120er3p0" "hltL1DoubleJet120er3p0" )  
# for runs 2017 C-v1,C-v2
if ( $era == "2017C-v1" || $era == "2017C-v2" ) then
   set ol1    = ( "hltL1DoubleJet40er3p0" "hltL1DoubleJet100er3p0" "hltL1DoubleJet112er3p0" "hltL1DoubleJet112er3p0" )
endif
set ocalo  = ( "hltDoubleCaloBJets30eta2p3" "hltDoubleCaloBJets100eta2p3" "hltDoubleCaloBJets100eta2p3" "hltDoubleCaloBJets100eta2p3" )
set obtag  = ( "hltBTagCalo30x8CSVp0p92SingleWithMatching" "hltBTagCalo80x6CSVp0p92SingleWithMatching" "hltBTagCalo80x6CSVp0p92SingleWithMatching" "hltBTagCalo80x6CSVp0p92SingleWithMatching" )
set opf    = ( "hltDoublePFJets40Eta2p3" "hltDoublePFJets100Eta2p3" "hltDoublePFJets200Eta2p3" "hltDoublePFJets350Eta2p3" )

@ count = 1
while ( $count < $#ptmin + 1 )
   set cfgfile = "tep_"$era"_"$ptmin[$count]"to"$ptmax[$count]".cfg"
   set outfile = "histograms_"$era"_tep_"$ptmin[$count]"to"$ptmax[$count]".root"
  
### CONFIG FILE

   echo $cfgfile

   cp -pf tag_and_probe.cfg $cfgfile
   
   sed -i -e "s/ERA/$era/g" $cfgfile
   sed -i -e "s/PTMIN/$ptmin[$count]/g" $cfgfile
   sed -i -e "s/PTMAX/$ptmax[$count]/g" $cfgfile
   sed -i -e "s/HLT/$hlt[$count]/g" $cfgfile
   sed -i -e "s/OL1/$ol1[$count]/g" $cfgfile
   sed -i -e "s/OCALO/$ocalo[$count]/g" $cfgfile
   sed -i -e "s/OBTAG/$obtag[$count]/g" $cfgfile
   sed -i -e "s/OPF/$opf[$count]/g" $cfgfile
   sed -i -e "s/JSON/$json/g" $cfgfile
   
   # btagcsv_2017C_ntuples.txt
   naf_submit.csh "btagcsv_"$era"_ntuples.txt" TagAndProbeOnlineBtag $cfgfile $splits $json
   
   rm -f $cfgfile

   @ count ++
end
