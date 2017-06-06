#!/bin/csh -f 


set sample = $1

cd $sample

set rootfiles = `find jobs/ -name "mssmhbb_triggers_*_x???.root"`

hadd "mssmhbb_triggers_"$sample"_control_back_final.root" $rootfiles


