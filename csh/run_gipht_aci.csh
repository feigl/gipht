#!/bin/csh -ex

# Run GIPhT from ACI
# 20140125 Kurt Feigl

#SBATCH --partition=geoscience          # queue 
#SBATCH --time=0-05:00:00		# run time in days-hh:mm:ss
#SBATCH --ntasks=16			# require 32 CPUs (CPUs)
#SBATCH --mem-per-cpu=4000		# RAM in MB (default 4GB, max 8GB)
#SBATCH --nodes=1                       # number of nodes requested
#SBATCH --ntasks-per-node=16            # default 16 if this line not specified
#SBATCH --cpus-per-task=1               # default 1 if this line not specified
#SBATCH --mem=32000                     # total RAM in MB, max 64GB  per node
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out


module purge
module list
module add comsol44
module add matlab-r2013b
module list
which matlab


/home/feigl/comsol_trials3/startcss

# this defaults to R2011b
#aci-005% ls -l /usr/local/bin/matlab
#lrwxrwxrwx 1 root root 48 Aug 28 11:35 /usr/local/bin/matlab -> /var/lib/condor/execute/MATLAB/R2011b/bin/matlab

#/usr/local/bin/matlab -nodisplay  <<EOF

# need release 2013b for comsol
#/usr/local/MATLAB/R2013b/bin/matlab -nodisplay  <<EOF
matlab -nodisplay  <<EOF
giphtpath
validate_comsol
%gipht
exit
EOF


/home/feigl/comsol_trials3/killcss

