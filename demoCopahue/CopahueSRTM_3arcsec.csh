#!/bin/csh -v

# for Copahue in Chile
 
### panels are labeled with Lower Left (SW) corner?
#foreach lat (N13 N14 N15)
#   foreach lon (W087 W088 W089 W090)

#panels are labeled with Lower Left (SW) corner?
foreach lat (S39 S38)
    foreach lon (W072 W071)
       echo Working on $lat $lon
       #ncftpget ftp://e0srp01u.ecs.nasa.gov/srtm/version2/North_America/$lat$lon.hgt.zip
       #wget ftp://e0srp01u.ecs.nasa.gov/srtm/version2/SRTM3/North_America/$lat$lon.hgt.zip
       wget -nd -v http://dds.cr.usgs.gov/srtm/version2_1/SRTM3/South_America/$lat$lon.hgt.zip
       unzip -v -o $lat$lon.hgt.zip
       srtm2grd.gmt $lat$lon.hgt
       echo done
   end # for lon
end # for lat

grdpaste -V S38W072.grd S38W071.grd -GS38row.grd
grdpaste -V S39W072.grd S39W071.grd -GS39row.grd
grdpaste -V S38row.grd  S39row.grd -GCopahueSRTM_3arcsec.grd
grdinfo CopahueSRTM_3arcsec.grd



# -T option means first line is at north edge
grd2xyz CopahueSRTM_3arcsec.grd -ZhT -V >! CopahueSRTM_3arcsec.i2
grdinfo CopahueSRTM_3arcsec.grd

grdsample CopahueSRTM_3arcsec.grd -R288.471665346/289.326385175/-38.237490935/-37.562215048 -I0.000277777/0.000277777 -GCopahueCSKdesc.grd
grdinfo CopahueCSKdesc.grd
grd2xyz CopahueCSKdesc.grd -ZhT -V >!  CopahueCSKdesc.i2 


