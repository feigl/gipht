#!/bin/csh 

#if ($#argv < 1) then
#cat - << ENDOFDOC
## Run GIPhT once for each pair listed in the list of interferograms specified 
#by ilist keyword in gipht.in
#ENDOFDOC
#
#else

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

# Now set up to run pairs individually
\rm -fv pair*.lst
set ilist = `grep ilist gipht.in | awk '$1 !~ /%/{print $3}'`
#echo $ilist
\cp gipht.in gipht.in.save

# make several files, each listing 1 interferometric pair
\rm -fv pair???.lst
# 2012-JUL-05 updated splitp.awk to start counting at pair001.lst
cat $ilist | grep a | awk -f `which splitp.awk`
#cat $ilist | grep b | awk -f `which splitp.awk`

foreach pair (pair???.lst)
echo '---------------------'
cat $pair

set imast = `cat $pair | awk 'NR==1{print $7}'`
set islav = `cat $pair | awk 'NR==1{print $8}'`

echo $pair $imast $islav

grep -v ilist gipht.in.save >! gipht.in
echo "ilist = $pair" >> gipht.in


matlab -nodisplay >! p_${imast}_${islav}.log  <<!
giphtpath
gipht
!

# name of output directory
set xdir = `ls -dt x_* | head -1`

\mv -v $xdir p_${imast}_${islav}_${xdir}

end # loop over pairs

\cp gipht.in.save gipht.in

#endif
