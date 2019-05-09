#!/bin/csh -f

if ( $#argv < 1 ) then
   echo "Need to give directory containing NAF jobs directories"
   exit
endif

set mydir = $1
cd $mydir
set mylist = `ls -1 | grep NAF_TagAndProbeOnlineBtag_tep`
if ( ! -d trash ) then
   mkdir trash
endif

foreach entry ( $mylist )
   set type = `echo $entry | awk -F "tep_" '{print $2}'`
   set hfile = "histograms_tep_"$type".root"
   if ( -d $entry) then
      if ( -e $hfile ) then
         rm -f $hfile
      endif
      hadd $hfile `find $entry -name "$hfile"` >>& log
      mv $entry trash
   endif
end

# # merge pt bins
# set eras = ( "QCD_Pt_50to80" "QCD_Pt_80to120" "QCD_Pt_120to170" "QCD_Pt_170to300" "QCD_Pt_300to470" "QCD_Pt_470to600" "QCD_Pt_600to800" "QCD_Pt_800to1000" "QCD_Pt_1000to1400" )
# set hfile = "histograms_tep_QCD.root"
# if ( -e $hfile ) then
#    rm -f $hfile
# endif
# foreach era ( $eras )
#    set myhfiles = `ls -1 "histograms_tep_"$era"_"*".root"`
#    hadd $hfile $myhfiles >>& log
# end

# # merge era versions
# if ( -e histograms_tep_2017C.root ) then
#    rm -f histograms_tep_2017C.root
# endif
# hadd histograms_tep_2017C.root histograms_tep_2017C-v1.root histograms_tep_2017C-v2.root histograms_tep_2017C-v3.root  >>& log
# if ( -e histograms_tep_2017E.root ) then
#    rm -f histograms_tep_2017E.root
# endif
# hadd histograms_tep_2017E.root  histograms_tep_2017E-v1.root histograms_tep_2017E-v2.root  >>& log
# 
# rm -f histograms_tep_2017CD.root histograms_tep_2017CDE-v1.root histograms_tep_2017CDE.root histograms_tep_2017E-v2F.root histograms_tep_2017CDEF.root >& /dev/null
# 
# hadd histograms_tep_2017CD.root histograms_tep_2017C.root histograms_tep_2017D.root  >>& log
# hadd histograms_tep_2017CDE-v1.root histograms_tep_2017CD.root histograms_tep_2017E-v1.root  >>& log
# hadd histograms_tep_2017CDE.root histograms_tep_2017CD.root histograms_tep_2017E.root  >>& log
# hadd histograms_tep_2017E-v2F.root histograms_tep_2017E-v2.root histograms_tep_2017F.root  >>& log
# hadd histograms_tep_2017CDEF.root histograms_tep_2017CDE.root histograms_tep_2017F.root  >>& log

sleep 1

touch trash
touch log
