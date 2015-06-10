#!/bin/csh -ex

# Run GIPhT 
# 20150609 Kurt Feigl

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
module add  matlab-r2014b
module list
which matlab

matlab -nodisplay  <<EOF
run(strcat(getenv('HOME'),filesep,'gipht',filesep,'src',filesep','giphtpath'));
%giphtpath
gipht
exit
EOF



