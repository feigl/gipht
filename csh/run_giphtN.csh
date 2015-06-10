#!/bin/csh 

#if ($#argv < 1) then
#cat - << ENDOFDOC
## Run GIPhT once for each pair listed in the list of interferograms specified 
#by ilist keyword in gipht.in
#ENDOFDOC
#
#else

echo Reminder 
echo 'RM -rf x_*'

\rm -fv pair*.lst


set ilist = `grep ilist gipht.in | awk '$1 !~ /%/{print $3}'`

echo $ilist

\cp gipht.in gipht.in.save


# make several files, each listing 1 interferometric pair
#grep a $ilist >! tmp.ilist 
#split -p a tmp.ilist ilist.
cat $ilist | grep a | awk -f `which splitp.awk`

foreach pair (pair*.lst)
echo '---------------------'
cat $pair

set imast = `cat $pair | awk 'NR==1{print $7}'`
set islav = `cat $pair | awk 'NR==1{print $8}'`

echo $pair $imast $islav

grep -v ilist gipht.in.save >! gipht.in
echo "ilist = $pair" >> gipht.in


matlab -nodisplay >! p${imast}_${islav}.log  <<!
giphtpath
gipht
!

# name of output directory
set xdir = `ls -dt x_* | head -1`

\mv -v $xdir p${imast}_${islav}_${xdir}

end

\cp gipht.in.save gipht.in

#endif
