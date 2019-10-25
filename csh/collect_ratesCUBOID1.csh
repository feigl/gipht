#!/bin/csh -vx
# Collect rates from output files from GIPHT produced by run_giphtN.csh
# 20190226 Kurt Feigl

cd CUBOID1
set out = 'dataCUBOID1.txt'

# get decimal years
grep -i time_fn_@_epoch_001_in_years____ p2*/x*.log | awk 'NF==9 && $4 == $9 {print $4}' >! tmp.epoch1
grep -i time_fn_@_epoch_002_in_years____ p2*/x*.log | awk 'NF==9 && $4 == $9 {print $4}' >! tmp.epoch2
#paste tmp.epoch1 tmp.epoch2 | awk '$1 == $3{print $1,$2}' >! tmp.epochs12
#paste tmp.epoch1 tmp.epoch2 | awk '{print $1,$2}' >! tmp.epochs12

# get YYYYMMDD
#grep In2 p2*/x*.log | grep 'Species A Member  2 Pair' | awk '{print $11}' | sed 's%/% %g' | awk '{print substr($9,3,8), substr($9,12,8)}' >! tmp.epochs12
grep time_fn_@_epoch_002_in_years____ p2*.log | sed 's%/% %g' | sed 's/p2/2/' | sed 's/_/ /' | sed 's/.log/ /' | awk '{print $1,$2}' | sort -u >! tmp.epochs12

# get estimated parameters
grep -i Volume_Change_in_cubic_meters___ p2*/x*.log | grep -v parameter | awk 'NF==3{print $0}' >! tmp.param_sigma
#grep -i SunDisk1_Excess_Pressure_in_Pa__ p2*/x*.log | grep -v parameter | sort | awk 'NF==3{print $0}' >! tmp.param_sigma

# write 1-line header
echo "Master  Slave Volume_Change_in_cubic_meters___ Uncertainty " >! $out               
#echo "Master  Slave SunDisk1_Excess_Pressure_in_Pa__ Uncertainty " >! $out               
# write data
#paste tmp.epochs12 tmp.param_sigma | sed 's/:F#//' | sed 's/:/ /' | awk '$1 == $4{print $2,$3,$6,$7}' >> $out               
#paste tmp.epochs12 tmp.param_sigma | sed 's/:F#//' | sed 's/:/ /' | awk '{print $1,$2,$5,$6}' >> $out                
#paste tmp.epoch1 tmp.epoch2 tmp.param_sigma | sed 's/:F#//' | sed 's/:/ /' | awk '{print $1,$2,$5,$6}' >> $out                
paste tmp.epochs12 tmp.param_sigma | awk '{print $1,$2,$4,$5}' >> $out                

