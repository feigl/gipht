#!/bin/csh
## Run GIPhT once for each pair listed in the list of interferograms specified by ilist keyword in gipht.in
# from ACI
# 20140716 Kurt Feigl
#SBATCH --partition=geoscience          # queue 
#SBATCH --time=0-24:00:00		# run time in days-hh:mm:ss
#SBATCH --ntasks=16			# require 32 CPUs (CPUs)
#SBATCH --mem-per-cpu=4000		# RAM in MB (default 4GB, max 8GB)
#SBATCH --nodes=1                       # number of nodes requested
#SBATCH --ntasks-per-node=16            # default 16 if this line not specified
#SBATCH --cpus-per-task=1               # default 1 if this line not specified
#SBATCH --mem=32000                     # total RAM in MB, max 64GB  per node
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

module add matlab-r2013b

/home/feigl/comsol_trials3/startcss

# save existing GIPHT file
set M = 0
while (! -e gipht.in.$M) 
  @ M = $M + 1
  \cp -v gipht.in gipht.in.$M
end

echo Reminder 
echo 'RM -rf x_*'
#\rm -rfv x_*

# first run ensemble
#matlab -nodisplay >! ENSEMBLE.log  <<!
matlab -nodisplay  << EOFE | tee  ENSEMBLE.log
giphtpath
addpath /usr/local/comsol44/mli -BEGIN
mphstart(2036)
gipht
EOFE


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
#cat $ilist | grep a | awk -f `which splitp.awk`
#cat $ilist | grep b | awk -f `which splitp.awk`
#on ACI, we need a local copy of this file. We have to take it with us.
cat $ilist | grep a | awk -f splitp.awk

foreach pair (pair???.lst)
echo '---------------------'
cat $pair

set imast = `cat $pair | awk 'NR==1{print $7}'`
set islav = `cat $pair | awk 'NR==1{print $8}'`

echo $pair $imast $islav

grep -v ilist gipht.in.save >! gipht.in
echo "ilist = $pair" >> gipht.in


#matlab -nodisplay >! p_${imast}_${islav}.log  <<!
matlab -nodisplay << EOF1 | tee p_${imast}_${islav}.log
giphtpath
addpath /usr/local/comsol44/mli -BEGIN
mphstart(2036)
gipht
EOF1

# name of output directory
set xdir = `ls -dt x_* | head -1`

\mv -v $xdir p_${imast}_${islav}_${xdir}

end # loop over pairs

\cp gipht.in.save gipht.in

#endif
#/home/feigl/comsol_trials3/killcss
