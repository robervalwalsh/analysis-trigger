#!/bin/csh -f

if ( $#argv < 7 ) then
   echo Need to give, config template, era, dataset, number of splits, trigger, ref trigger, json
   exit
endif

set cfgin = $1
set era = $2
set dataset = $3
set splits = $4
set trigger = $5
set triggerref = $6
set json = $7

   set cfgbase = "pfjet_"$dataset"_"$era"_"$trigger"_ref"$triggerref
   set cfgfile = $cfgbase".cfg"
   set outfile = "histograms_"$cfgbase".root"
  
### CONFIG FILE

   echo $cfgfile

   cp -pf $cfgin $cfgfile
   
   sed -i -e "s/ERA/$era/g" $cfgfile
   sed -i -e "s/OUTFILE/$outfile/g" $cfgfile
   sed -i -e "s/JSON/$json/g" $cfgfile
   
   naf_submit.csh $dataset"_"$era"_ntuples.txt" SingleJetTriggerEfficiencies $cfgfile $splits $json
   
   rm -f $cfgfile
