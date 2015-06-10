#!/bin/csh 

# Run GIPhT from ACI
# 20140125 Kurt Feigl

set HERE = $PWD

#foreach ING (T344_0.ing)
#foreach ING (*.ing)
#foreach ING (ENV_T115.ing)
foreach ING (ENV_T115.ing  ENV_T344.ing  ERS_T072.ing  ERS_T115.ing  ERS_T179.ing  ERS_T344.ing  ERS_T451.ing)
#echo '---------------------'
#echo '---------------------'
#echo $ING
#echo '---------------------'
set TAG = $ING:r

set DATAFILENAME = `grep datafilename $ING | grep -v % | awk '{print $3}'`
#echo DATAFILENAME is $DATAFILENAME

if (! -e ../gipht28.$TAG) then
  mkdir ../gipht28.$TAG
endif

cd ../gipht28.$TAG
\cp -f $HERE/$TAG.ing gipht.in
\rm -f *.gin          ;  ln -s ../gipht28/*.gin .
\rm -f giphtpath.m    ;  ln -s ../gipht28/giphtpath.m .
\rm -f *.awk          ;  ln -s ../gipht28/*.awk .
\rm -f *.csh          ;  ln -s ../gipht28/*.csh .

# get data file for comsol
if ($DATAFILENAME != "") then
if (-e $HERE/$DATAFILENAME ) then
   \cp -f $HERE/$DATAFILENAME .
endif
endif


#echo sbatch -J $TAG -p geoscience run_gipht_aci.csh
echo cd ../gipht28.$TAG  
echo sbatch -J $TAG -p geoscience run_gipht_n_pairs_aci.csh 
#echo cd $HERE


cd $HERE
end # loop
#

