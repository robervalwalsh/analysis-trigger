#!/bin/csh -f

if ( $#argv < 5 ) then
   echo Need to give, config template, era, dataset, number of splits, trigger
   exit
endif

set cfgin = $1
set era = $2
set dataset = $3
set splits = $4
set trigger = $5

#set ptmin =  (20  100 200 300 350 650 )
#set ptmax =  (100 200 300 450 650 1000 )

set ptmin =  (60  )
set ptmax =  (500 )
set etamax = 2.2

set refpt = (40 60 80 140 200 260)

# for HLT_PFJet100
if ( $trigger == "HLT_PFJet100" ) then
   set l1ptmin = 100
   set l1etamax = 2.3
   set l2ptmin = 100
   set l2etamax = 2.3
   set l3ptmin = 100
   set l3etamax = 2.3
else if ( $trigger == "HLT_PFJet40" ) then
   set l1ptmin = 0
   set l1etamax = 3.0
   set l2ptmin = 10
   set l2etamax = 3.0
   set l3ptmin = 40
   set l3etamax = 3.0
else if ( $trigger == "HLT_PFJet60" ) then
   set l1ptmin = 35
   set l1etamax = 3.0
   set l2ptmin = 40
   set l2etamax = 3.0
   set l3ptmin = 60
   set l3etamax = 3.0
else if ( $trigger == "HLT_PFJet80" ) then
   set l1ptmin = 60
   set l1etamax = 3.0
   set l2ptmin = 50
   set l2etamax = 3.0
   set l3ptmin = 80
   set l3etamax = 3.0
else if ( $trigger == "HLT_PFJet140" ) then
   set l1ptmin = 90
   set l1etamax = 3.0
   set l2ptmin = 110
   set l2etamax = 3.0
   set l3ptmin = 140
   set l3etamax = 3.0
else if ( $trigger == "HLT_PFJet200" ) then
   set l1ptmin = 120
   set l1etamax = 3.0
   set l2ptmin = 170
   set l2etamax = 3.0
   set l3ptmin = 200
   set l3etamax = 3.0
else if ( $trigger == "HLT_PFJet260" ) then
   set l1ptmin = 200
   set l1etamax = 3.0
   set l2ptmin = 210
   set l2etamax = 3.0
   set l3ptmin = 260
   set l3etamax = 3.0
endif 



set hlt    = ( "HLT_PFJet60_v"      "HLT_PFJet80_v"      "HLT_PFJet140_v"      "HLT_PFJet200_v"      "HLT_PFJet260_v"  )
set l1     = ( "L1_SingleJet35"     "L1_SingleJet60"     "L1_SingleJet90"      "L1_SingleJet120"     "L1_SingleJet200" )
set ol1    = ( "hltL1sSingleJet35"  "hltL1sSingleJet60"  "hltL1sSingleJet90"   "hltL1sSingleJet120"  "hltL1sSingleJet170IorSingleJet180IorSingleJet200" )
set ocalo  = ( "hltSingleCaloJet40" "hltSingleCaloJet50" "hltSingleCaloJet110" "hltSingleCaloJet170" "hltSingleCaloJet210" )
set opf    = ( "hltSinglePFJet60"   "hltSinglePFJet80"   "hltSinglePFJet140"   "hltSinglePFJet200"   "hltSinglePFJet260" )

@ count = 1
while ( $count < $#ptmin + 1 )
   set cfgbase = "pfjet_"$dataset"_"$era"_hlt"$l3ptmin"_refhlt"$refpt[$count]"_pt"$ptmin[$count]"to"$ptmax[$count]
   set cfgfile = $cfgbase".cfg"
   set outfile = "histograms_"$cfgbase".root"
  
### CONFIG FILE

   echo $cfgfile

   cp -pf $cfgin $cfgfile
   
   sed -i -e "s/ERA/$era/g" $cfgfile
   sed -i -e "s/OUTFILE/$outfile/g" $cfgfile
   sed -i -e "s/OFFPTMIN/$ptmin[$count]/g" $cfgfile
   sed -i -e "s/OFFPTMAX/$ptmax[$count]/g" $cfgfile
   sed -i -e "s/OFFETAMAX/$etamax/g" $cfgfile
   sed -i -e "s/MyHLT/$hlt[$count]/g" $cfgfile
   sed -i -e "s/MyL1/$l1[$count]/g" $cfgfile
   sed -i -e "s/OL1/$ol1[$count]/g" $cfgfile
   sed -i -e "s/OCALO/$ocalo[$count]/g" $cfgfile
   sed -i -e "s/OPF/$opf[$count]/g" $cfgfile
   sed -i -e "s/L1PTMIN/$l1ptmin/g" $cfgfile
   sed -i -e "s/L1ETAMAX/$l1etamax/g" $cfgfile
   sed -i -e "s/L2PTMIN/$l2ptmin/g" $cfgfile
   sed -i -e "s/L2ETAMAX/$l2etamax/g" $cfgfile
   sed -i -e "s/L3PTMIN/$l3ptmin/g" $cfgfile
   sed -i -e "s/L3ETAMAX/$l3etamax/g" $cfgfile
   
   naf_submit.csh $dataset"_"$era"_ntuples.txt" SingleJetTriggerAnalysis $cfgfile $splits json_2017.txt
   
   rm -f $cfgfile

   @ count ++
end
