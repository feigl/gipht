function [emin,emax,nmin,nmax] = select_subset(siteCode,subset)
%function [emin,emax,nmin,nmax] = select_subset(siteCode,subset)
% given siteCode and geographic subset, 
%return bounding box in UTM coordinates in meters
% 2021/06/22 Kurt Feigl
switch siteCode
    case 'COSOC'
        switch subset
            case 'all'
                emin = -Inf;
                emax = +Inf;
                nmin = -Inf;
                nmax = +Inf;
            case 'big'
                emin =  418000;
                emax =  436000;
                nmin = 3980000;
                nmax = 3998000;
            case 'med'  % square
                emin =  418000;
                emax =  434000;
                nmin = 3980000;
                nmax = 3996000;
            case 'sml'  %
                emin = 294500;
                emax = 297500;
                nmin = 4471000;
                nmax = 4475000;
            otherwise
                error(sprintf('unrecognized subset %s\n',subset));
        end
    case 'SANEM'
        switch subset
            case 'all'
                emin = -Inf;
                emax = +Inf;
                nmin = -Inf;
                nmax = +Inf;
           case 'big'
                emin = nanmin(TMAP.Easting);
                emax = nanmax(TMAP.Easting);
                nmin = nanmin(TMAP.Northing);
                nmax = nanmax(TMAP.Northing);
            case 'med'  % Coincide with maps made by Matt Folsom
                emin = 293500;
                emax = 301500;
                nmin = 4466500;
                nmax = 4478500;
            case 'sml'  %
                emin = 294500;
                emax = 297500;
                nmin = 4471000;
                nmax = 4475000;
            otherwise
                error(sprintf('unrecognized subset %s\n',subset));
        end
    otherwise
        error(sprintf('Unknown siteCode %s\n',siteCode'));       
end
return
end

