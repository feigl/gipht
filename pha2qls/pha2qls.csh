#!/bin/csh


# run pha2qls.c to perform quad-tree resampling on a phase file
if ($#argv < 1) then
cat - << ENDOFDOC
$0:t Last modification 20150819 Kurt Feigl

perform Quad-tree resampling on wrapped phase
usage pha2qls [options] input.pha nx ny 
Required arguments: 
        phase.grd == Name of input binary file containing phase values coded as one signed byte per such that 256 DN = 1 cycle and values range from -128 to +127
Output: 
        File containing Quadtree List of Samples (QLS) 6-column binary file containing I,J,W,Q,GX,GY [named output.qls by default]
Options: 
        -h   this help information 
        -L l == threshold Level   for circular mean deviation misfit to 1-parameter mean [256 DN = 1 cycle of phase] 
        -M m == threshold Maximum for circular mean deviation misfit to 3-parameter ramp [256 DN = 1 cycle of phase] 
        -N n == mininum number of pixels in patch 
        -o output_qls_filename.pha   == file name for output reconstructed phase values, written as an image with nx columns and ny rows 
        -P output_phase_filename.pha == file name for output reconstructed phase values, written as an image with nx columns and ny rows.
        -X output_gradx_filename.i2  == file name for differential phase in positive X-direction (increasing column index) [2^16 = 65536 DN = 1 cycle of phase per pixel]
        -Y output_gradx_filename.i2  == file name for differential phase in negative Y-direction (increasing line index)   [2^16 = 65536 DN = 1 cycle of phase per pixel]
        -V Verbose mode
        -d n  optional debug levels
                 -d 0  == No additional output debugging information to stderr (default)
                 -d 1  == Warning to stdout if a patch write exceeds nx,ny bounds or if a patch write overwrites an existing patch
                 -d 2  == echos input qls data to stdout
                 -d 3  == combines 1 & 2

ENDOFDOC

exit 1
endif

# initialiaze
set CURRENTDIR = `\pwd`
set Nminpix = 16         
set Mmaxcmd = 8  
set Levelmx = 127
set isll    = 1    # == 1 if geographic coordinates (longitude, latitude in degrees)
set rapixm = 1.;  # dimension of range pixel in meters   - do not know this yet
set azpixm = 1.;  # dimension of azimuth pixel in meters - do not know this yet
set isverbose = 0
set VERBOSE = ''

#


# parse command line
while ($#argv > 0 )
   set input = ( $argv )
   switch($input[1])
        case -i:
                set INFILE  = $input[2]
                shift input
                breaksw
        case -o:
                set qhase  = $input[2]
                shift input
                breaksw
        case -L:
                set Levelmx  = $input[2]
                shift input
                breaksw
        case -M:
                set Mmaxcmd  = $input[2]
                shift input
                breaksw
        case -N:
                set Nminpix  = $input[2]
                shift input
                breaksw
        case -R:
                set isll     = 0;
                shift input
                breaksw
        case -V:
                set isverbose     = 1;
                set VERBOSE = '-V'
                shift input
                breaksw

   #      case -p:  
   #              set DOp = 1           
   #              set PICKPAIROPT1 = ( ${input[2-]} )
   #              set PICKPAIROPT  = `echo $PICKPAIROPT1 |awk '{printf("-O @ %s @ \n",$0)}' | sed "s/@/\'/g"`
   #              echo Setting PICKPAIROPT to $PICKPAIROPT  
   #              breaksw
    endsw
   shift argv
end

# document
set maker = $0:t
set today = `date`
set remarked = `echo by $USER on $today with $maker`
if ($isverbose) echo remarked is $remarked

# find number of rows
set nrows = `gmt grdinfo $INFILE -C | awk '{print $11}'`
if ($isverbose) echo nrows = $nrows

# find number of columns
set ncols = `gmt grdinfo $INFILE -C | awk '{print $10}'`
if ($isverbose) echo ncols = $ncols

# find ranges
set ranges = `gmt grdinfo $INFILE -Ir | head -1`
if ($isverbose) echo ranges is $ranges

# find increments
set incx = `gmt grdinfo $INFILE -Vq | grep x_inc | awk '{print $7}'`
if ($?incx == 0) then
  set incx = "1"
endif
if ($isverbose) echo incx is $incx
set incy = `gmt grdinfo $INFILE -Vq | grep y_inc | awk '{print $7}'`
if ($?incy == 0) then
  set incy = "1"
endif
if ($isverbose) echo incy is $incy

# find mid latitude
set midlat = `gmt grdinfo $INFILE -C | awk '{print ($4 + $5)/2.}'`
if ($isverbose) echo midlat is $midlat

# find mid longitude
set midlon = `gmt grdinfo $INFILE -C | awk '{print ($2 + $3)/2.}'`
if ($isverbose) echo midlon is $midlon

# find the name of the parameter file 
set prmfile = `find . -type f -name '*.PRM' -print | awk 'NR ==1 {print $1}' `
if ($isverbose) echo prmfile is $prmfile

# find fringe spacing in m
set metersperfringe = `grep radar_wavelength $prmfile | awk 'NR == 1 {print $3/2.}'`
if ($isverbose) echo metersperfringe is $metersperfringe

# try to decide if file is geographic (lat, lon) coordinates
gmt grdinfo $INFILE | grep -e Transverse_Mercator -e UTM >! tmp.utm
echo $INFILE | grep "_ll" >! tmp.ll

if (! -z tmp.ll && -z tmp.utm) then  # Geographic coordinates (lat, longitude) in degrees
  set icoord = 1
  set fswitch = "-fo0g,1g"
  set FACTORX = "XINC DEG2KM 1000 MUL"
  set FACTORY = "YINC DEG2KM 1000 MUL"
  set XNAME = "longitude in degrees_east"
  set YNAME = "latitude in degrees_north"
else 
  if (! -z tmp.utm) then # UTM
    set iscoord = 2
    set fswitch = "-fo0f,1f"
    set FACTORX = "XINC"
    set FACTORY = "YINC"
    set XNAME = "UTM Easting in meters"
    set YNAME = "UTM Northing in meters"
  else # pixel indices in range, azimuth
   set iscoord = 0
   set fswitch = "-fo0f,1f"
    set FACTORX = "1"
    set FACTORY = "1"
    set XNAME = "IX"
    set YNAME = "IY"
  endif
endif

if ($isverbose) echo INFILE is $INFILE
if ($isverbose) echo XNAME is $XNAME
if ($isverbose) echo YNAME is $YNAME
if ($isverbose) echo FACTORX is $FACTORX
if ($isverbose) echo FACTORY is $FACTORY

# convert from radians to DN
gmt grdmath $INFILE ISFINITE $INFILE MUL PI DIV 2.0 DIV 256 MUL = tmp.grd
if ($isverbose) gmt grdinfo -L2 tmp.grd
gmt grd2xyz tmp.grd -ZTLc >! tmp.pha

# make output file names by defaults
set qhase = $INFILE:r\_qphase.grd
set gradx = $INFILE:r\_qgradx.grd
set grady = $INFILE:r\_qgrady.grd
set gradt = $INFILE:r\_qgradt.grd
set qmask = $INFILE:r\_qmaskn.grd

# perform Quad-tree resampling on wrapped phase
# usage pha2qls [options] input.pha nx ny 
# Required arguments: 
#         input.pha == Name of input binary file containing phase values coded as one signed byte per such that 256 DN = 1 cycle and values range from -128 to +127
#         nx        == Number of pixels in X-direction [= number of columns]
#         ny        == Number of pixels in Y-direction [= number of rows]
# Output: 
#         File containing Quadtree List of Samples (QLS) 6-column binary file containing I,J,W,Q,GX,GY [named output.qls by default]
# Options: 
#         -h   this help information 
#         -L l == threshold Level   for circular mean deviation misfit to 1-parameter mean [256 DN = 1 cycle of phase] 
#         -M m == threshold Maximum for circular mean deviation misfit to 3-parameter ramp [256 DN = 1 cycle of phase] 
#         -N n == mininum number of pixels in patch 
#         -o output_qls_filename.pha   == file name for output reconstructed phase values, written as an image with nx columns and ny rows 
#         -P output_phase_filename.pha == file name for output reconstructed phase values, written as an image with nx columns and ny rows.
#         -X output_gradx_filename.i2  == file name for differential phase in positive X-direction (increasing column index) [2^16 = 65536 DN = 1 cycle of phase per pixel]
#         -Y output_gradx_filename.i2  == file name for differential phase in negative Y-direction (increasing line index)   [2^16 = 65536 DN = 1 cycle of phase per pixel]
#         -V Verbose mode
#         -d n  optional debug levels
#                  -d 0  == No additional output debugging information to stderr (default)
#                  -d 1  == Warning to stdout if a patch write exceeds nx,ny bounds or if a patch write overwrites an existing patch
#                  -d 2  == echos input qls data to stdout
#                  -d 3  == combines 1 & 2

if (! -e /Users/feigl/gipht/pha2qls/pha2qls.maci64 ) then
  echo "ERROR: cannot find executable for pha2qls.c"
  exit 1
endif


/Users/feigl/gipht/pha2qls/pha2qls.maci64 tmp.pha $ncols $nrows \
$VERBOSE -P quadt.pha -X gradx.i2 -Y grady.i2 -T qls.txt \
-L $Levelmx -N $Nminpix -M $Mmaxcmd -Q 10000 >&! pha2qls.log

# convert output files to grids
if ($isverbose) echo ranges is $ranges 
if ($isverbose) echo incx is $incx
if ($isverbose) echo incy is $incy

# -Dxname/yname/zname/scale/offset/invalid/title/remark 
# wrapped phase after quad-tree
gmt xyz2grd quadt.pha $ranges -I$incx/$incy -di0 -r -ZTLc $fswitch -Gtmp.grd 
gmt grdmath tmp.grd 256 DIV 2.0 MUL PI MUL = $qhase
gmt grdedit -D"$XNAME"/"$YNAME"/"radians"/1///"wrapped phase after quadtree"/"$remarked" $qhase
if ($isverbose) gmt grdinfo -L2 $qhase

# east component of gradient
gmt xyz2grd gradx.i2  $ranges -I$incx/$incy -di0 -r -ZTLh $fswitch -Gtmp.grd
gmt grdmath tmp.grd 256 DIV 256 DIV $metersperfringe MUL $FACTORX DIV = $gradx 
gmt grdedit -D"$XNAME"/"$YNAME"/"dimensionless"/1///"east component of gradient"/"$remarked" $gradx
if ($isverbose) gmt grdinfo -L2 $gradx

# north component of gradient
gmt xyz2grd grady.i2  $ranges -I$incx/$incy -di0 -r -ZTLh $fswitch -Gtmp.grd
gmt grdmath tmp.grd 256 DIV 256 DIV $metersperfringe MUL $FACTORY DIV = $grady 
gmt grdedit -D"$XNAME"/"$YNAME"/"dimensionless"/1///"north component of gradient"/"$remarked" $grady
if ($isverbose) gmt grdinfo -L2 $grady

# mask containing quadtree value at centroids of quadtree patches and NaN elsewhere
if ($isverbose) echo nrows = $nrows
if ($isverbose) echo ncols = $ncols
cat qls.txt | awk 'NR > 1{print $1,$2,$4}' >! tmp.ijq
#gmt gmtinfo tmp.ijq


#gmt nearneighbor tmp.ijq -S1 -N1 -R0/$ncols/0/$nrows -I1/1 -r -f0f,1f,2f | gmt xyz2grd -R0/$ncols/0/$nrows -I1/1 -r -Gijq.grd
gmt blockmedian tmp.ijq -C -R0/$ncols/0/$nrows -I1/1 -r -f0f,1f,2f | gmt xyz2grd -R0/$ncols/0/$nrows -I1/1 -r -Gijq.grd
if ($isverbose) gmt grdinfo -L2 ijq.grd

# normalize mask 
gmt grdmath ijq.grd ijq.grd DIV = ijm.grd 
if ($isverbose) gmt grdinfo -L2 ijm.grd

# copy values into coordinates of input gird
gmt grd2xyz ijm.grd -ZTLd | gmt xyz2grd -ZTLd $ranges -I$incx/$incy $fswitch -G$qmask -r
gmt grdedit -D"$XNAME"/"$YNAME"/"mask"/1///"quadtree mask"/"$remarked" $qmask
if ($isverbose) gmt grdinfo -L2 $qmask

\rm -f tmp.ijq ijq.grd ijm.grd 





