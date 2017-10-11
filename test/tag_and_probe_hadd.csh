#!/bin/csh -f

if ( $#argv < 1 ) then
   echo Need to give era
   exit
endif

set era = $1

set ptmin =  (40  100 200 350)
set ptmax =  (100 200 350 1000)

@ count = 1
set haddall = ""
while ( $count < $#ptmin + 1 )
   set nafdir = "NAF_TagAndProbeOnlineBtag_tep_"$era"_"$ptmin[$count]"to"$ptmax[$count]
   set hfile = "histograms_"$era"_tep_"$ptmin[$count]"to"$ptmax[$count]
   set hdirfile = $nafdir"/"$hfile"_all.root"
   set haddall = "$haddall $hdirfile"
   cd $nafdir
   if ( -e $hfile"_all.root" ) then
      rm -f $hfile"_all.root"
   endif
   hadd $hfile"_all.root" `find . -name "$hfile.root"`
   cd ..
   @ count ++
end

set haddfinal = "histograms_"$era"_tep_"$ptmin[1]"to"$ptmax[$#ptmax]"_all.root"
if (  -e $haddfinal ) then
   rm -f $haddfinal 
endif

hadd $haddfinal $haddall
