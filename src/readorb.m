function [xs, ys, zs, xdot, ydot, zdot, mjd1950, sec, orbnum] = readorb(orbfile)
% read a diapason format .orb file in m, m/s and s

% MJD1950 is apparently the number of days since 1950/01/01
%
% Date 1950/01/01  0:00 hrs, DOY   1 JD  2433282.5000 MJD  33282.0000
% GPS Week -1567 Day of week  0, GPS Seconds      0 Day of Week Sun
% Decimal Year  1950.00000
% hengill.geology.wisc.edu% doy 1985 1 1
% Date 1985/01/01  0:00 hrs, DOY   1 JD  2446066.5000 MJD  46066.0000
% GPS Week   260 Day of week  2, GPS Seconds 172800 Day of Week Tue
% Decimal Year  1985.00000
% hengill.geology.wisc.edu% doy 0
% Date 1858/11/17  0:00 hrs, DOY 321 JD  2400000.5000 MJD      0.0000
% GPS Week -6321 Day of week  3, GPS Seconds 259200 Day of Week Wed
% Decimal Year  1858.87671
% hengill.geology.wisc.edu% doy 48736
% Date 1992/04/24  0:00 hrs, DOY 115 JD  2448736.5000 MJD  48736.0000
% GPS Week   641 Day of week  5, GPS Seconds 432000 Day of Week Fri
% Decimal Year  1992.31148

if fexist(orbfile)    
    [mjd1950,sec,orbnum,xs,ys,zs,xdot,ydot,zdot] = textread(orbfile,'%d%f%d%*[^\n]\n%f%f%f\n%f%f%f\n',-1);
    
    n = max([numel(sec),numel(xs),numel(xdot)]);
    fprintf (1,'Read %5d records from orbit file %s\n',n,orbfile);
    
    sec = sec/1000.0; % convert from milliseconds to seconds
    xs  = xs*1000.0;  % convert from km to m
    ys  = ys*1000.0;  % convert from km to m
    zs  = zs*1000.0;  % convert from km to m
    xdot = xdot*1000.0; % convert from km/s to m/s
    ydot = ydot*1000.0; % convert from km/s to m/s
    zdot = zdot*1000.0; % convert from km/s to m/s
    
    if abs(max(mjd1950) - min(mjd1950)) > 0
        warning(sprintf('MJD1950 changes during pass %d %d %d\n',min(mjd1950),max(mjd1950),max(mjd1950)-min(mjd1950)));
    end
    if abs(max(orbnum) - min(orbnum)) > 0
        warning(sprintf('Orbit number changes during pass %d %d %d\n',min(orbnum),max(orbnum),max(orbnum) - min(orbnum)));
    end
    if min(sec) < 0
        warning(sprintf('Number of seconds is negative %#E12.6 in %s\n',min(sec),orbfile));
    end
    if max(sec) > 86400.0
        warning(sprintf('Number of seconds exceeds 86400 %#E12.6 in %s\n',min(sec),orbfile));
    end
    if min(mjd1950)  < 14610
        warning(sprintf('mjd1950 (%#12.6E) is less than 14610 (before 1990-JAN-01) in %s\n',min(mjd1950),orbfile));
    end
    if max(mjd1950) > 25567
        warning(sprintf('mjd1950 (%#12.6E) is greater than 25567 (after 2020-JAN-01) in %s\n',max(mjd1950),orbfile));
    end
    
    iok = find(sec>0.0);
    iok = intersect(iok,find(sec<=86400.0));
    iok = intersect(iok,find(mjd1950>14610));
    iok = intersect(iok,find(mjd1950<25567));
    if n - numel(iok)  > 0
        warning(sprintf('Pruning %d invalid records from %s\n',n - numel(iok),orbfile));
        mjd1950=mjd1950(iok);
        sec=sec(iok);
        orbnum=orbnum(iok);
        xs=xs(iok);
        ys=ys(iok);
        zs=zs(iok);
        xdot=xdot(iok);
        ydot=ydot(iok);
        zdot=zdot(iok);
    else
        fprintf (1,'Returning %5d valid records from orbit file %s\n',numel(iok),orbfile);
    end
else
    error(sprintf('Could not find orbit file named %s\n',orbfile));
end

return;
