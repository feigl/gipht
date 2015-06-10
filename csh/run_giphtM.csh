#!/bin/csh


clean
module add matlab-r2013b
#-rw-rw-r--  1 feigl feigl 6.9k Dec  8 19:19 gipht.in.just2
#-rw-rw-r--  1 feigl feigl 7.1k Dec  8 19:47 gipht.in.pgradc20081
#-rw-rw-r--  1 feigl feigl 7.0k Dec  8 19:53 gipht.in.Post1997
#-rw-rw-r--  1 feigl feigl 6.9k Dec  8 19:52 gipht.in.Pre1997
#-rw-rw-r--  1 feigl feigl 6.2k Dec  7 13:53 gipht.in.save
#
#foreach infile (gipht.in.*)
#foreach infile (gipht.in.pgradc20081 gipht.in.NoCo)
#foreach infile (gipht.in.pgradc1997 gipht.in.Post1997tf  )
#foreach infile (gipht.in.Pgradc1997 gipht.in.Pgradc2008)
#foreach infile (gipht.in.Co1997 gipht.in.Co2008)
#foreach infile (gipht.in.Mogi1997 gipht.in.Mogi2008 gipht.in.Pgradc1997 gipht.in.Pgradc2008)

## Run Ensemble ##

#foreach infile (*.ing)
#\cp -v $infile gipht.in

#matlab -nodisplay <<!
#giphtpath
#gipht
#!

#set rn = `ls -1dt x_* | awk 'NR==1{print $1}'`

#\mv -v $rn  $infile:r_$rn

#end # loop over files

# save existing GIPHT file
set M = 0
while (! -e gipht.in.$M)
  @ M = $M + 1
  \cp -v gipht.in gipht.in.$M
end

echo Reminder
echo 'RM -rf x_*'
\rm -rfv x_*

# first run ensemble
matlab -nodisplay >! ENSEMBLE.log  <<!
giphtpath
gipht
!

# name of output directory
set xdir = `ls -dt x_* | head -1`
\mv -v $xdir ENSEMBLE_${xdir}
\mv -v ENSEMBLE.log ENSEMBLE_${xdir}





