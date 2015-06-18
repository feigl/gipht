#!/bin/csh

# convert the file into 4-byte float
imageMath.py -e='a_1' --a=filt_topophase.unw.geo -o phase.geo -t FLOAT

# convert the file into 4-byte float
imageMath.py -e='a_0' --a=filt_topophase.unw.geo -o ampli.geo -t FLOAT

# convert the file into 4-byte float
imageMath.py -e='a_0' --a=phsig.cor.geo          -o coher.geo -t FLOAT

foreach part (phase ampli coher)
echo working on part $part

# make an envi header
isce2gis.py envi -i $part.geo
  
# make an vert header
isce2gis.py vrt -i $part.geo

# get number of columns
set ncols = `grep samples $part.geo.hdr | awk '{print $3}'`
echo number of colums is $ncols

# get number of rows
set nrows = `grep lines $part.geo.hdr | awk '{print $3}'`
echo number of rows is $nrows

# get lon1
set lon1 = `grep GeoTransform $part.geo.vrt | sed 's/>/,/g' | sed 's/</,/g' | awk -F, '{printf("%.10f\n", $3)}'`
echo lon1 is $lon1

# get dlon
set dlon = `grep GeoTransform $part.geo.vrt | sed 's/>/,/g' | sed 's/</,/g' | awk -F, '{printf("%.10f\n",  $4)}'`
echo dlon is $dlon

# get abs(dlon)
set absdlon = `echo $dlon | awk -F, '{printf("%.10f\n",sqrt($1*$1))}'`
echo dlat is $absdlon

# calculate lon2
set lon2 = `echo $lon1 $dlon $ncols | awk '{printf("%.10f\n",$1+$2*($3-1))}'`
echo lon2 is $lon2

# calculate lonW and lonE
set lonWE = `echo $lon1 $lon2 | awk '{if ($1 < $2) {printf("%.10f/%.10f\n",$1,$2)} else {printf("%16.10f/%16.10f\n",$2,$1)}}'`
echo lonWE is $lonWE

# get lat1
set lat1 = `grep GeoTransform $part.geo.vrt | sed 's/>/,/g' | sed 's/</,/g' | awk -F, '{printf("%.10f\n",  $6)}'`
echo  lat1 is $lat1

# get dlat
set dlat = `grep GeoTransform $part.geo.vrt | sed 's/>/,/g' | sed 's/</,/g' | awk -F, '{printf("%.10f\n",  $8)}'`
echo dlat is $dlat

# get abs(dlat)
set absdlat = `echo $dlat | awk -F, '{printf("%.10f\n", sqrt($1*$1))}'`
echo dlat is $absdlat

# calculate lat2
set lat2 = `echo $lat1 $dlat $nrows | awk '{printf("%.10f\n",$1+$2*($3-1))}'`
echo lat2 is $lat2

# calculate latS and latN
set latSN = `echo $lat1 $lat2 | awk '{if ($1 < $2) {printf("%.10f/%.10f\n",$1,$2)} else {printf("%.10f/%.10f\n",$2,$1)}}'`
echo latSN is $latSN

# Read 4-byte floating point single precision from Left and Top
xyz2grd -V -R$lonWE/$latSN -I$absdlon/$absdlat -ZTLf $part.geo -G$part.grd

# check
grdinfo $part.grd

end # loop over parts







