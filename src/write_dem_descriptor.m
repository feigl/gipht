function ierr = write_dem_descriptor(demdescfile,demdescfile2,isgeo,y1,x1,nl,nc,l1,c1,ml,mc,dl,dc,fi2,lat0,lon0,y0,x0,hemisphere,iutmzone,isext)
%function ierr = write_dem_descriptor(demdescfile,demdescfile2,isgeo,y1,x1,nl,nc,l1,c1,ml,mc,dl,dc,fi2,lat0,lon0,y0,x0,hemisphere,iutmzone,isext)
% write DEM parameters to descriptor file
% input:
%     demdesc  == file name of DEM descriptor, e.g. 'dem.dat' points to 'dem.i2'
%     demdesc 2 == file name of DEM descriptor, e.g. 'dem.dat' points to 'dem.i2'
%      isgeo   == 1 for geographic coordinates in degrees
%              == 2 for UTM coordinates in meters
%              == 3 for Lambert coordinates in meters
%      y1,x1   == line (latitude), column (longitude) coordinates of first pixel
%      nl, nc  == total numbers of lines and columns
%      l1, c2  == line and col coordinates of first pixel to extract
%      ml, mc  == number of columns to extract
%      dl, dc  == step sizes in lines and columns [MAY BE NEGATIVE]
%      fi2     == file name of binary I2 DEM file e.g., 'dem.i2'
%
%for use with DIAPASON DTOOLS
%
% Modifications
% 2009-MAR-30 Kurt admit blank lines
% 2009-APR-01 Do some more checking
% 2009-JUL-20 Return cartographic origin lon0, lat0
% 2010-APR-02 Handle Lambert

fprintf(1,'%s begins ...\n',mfilename);

ierr = 100;
isext = 1;

fid = fopen(demdescfile,'r');
if fid <= 0
    fprintf (1,'ERROR: Cannot open DEM descriptor file called %s\n',demdescfile);
    ierr = 1;
    error 'Cannot open DEM descriptor'
    return
end
fid2 = fopen(demdescfile2,'w+');
if fid2 <= 0
    fprintf (1,'ERROR: Cannot open DEM descriptor file called %s\n',demdescfile2);
    error 'Cannot open DEM descriptor'
    ierr = 2;
    return
end

%line1 = 0;
while 1 % for line in the file
    % read one line and make sure it is a string
    line1 = fgetl(fid);
    %line1 = fscanf(fid,'%s\n')
    %if ~isstr(line1)
    %if line1 == -1
    if ~ischar(line1)
        disp 'Reached end of DEM descriptor file'
        fclose(fid);
        fclose(fid2);
        %         if (isext == 0)
        %             ml = nl;
        %             mc = nc;
        %             l1 = 1;
        %             c1 = 1;
        %         else
        %             if l1 + ml > nl ||  c1 + mc > nc ||  l1 < 1 ||  c1 < 1 ||  mc < 1 ||  ml < 1 ||  nl < 1 ||  nc < 1
        %                 l1
        %                 ml
        %                 nl
        %                 c1
        %                 mc
        %                 nc
        %                 error ('ERROR: invalid dimensions for DEM');
        %             end
        %         end
        ierr = 0;
        break
    else
        if numel(strfind(line1,'EXTRACTION')) > 0 && numel(strfind(line1,'PAS')) == 0
            if isext == 1
                fprintf(fid2,'EXTRACTION %s\n','OUI');
            else
                fprintf(fid2,'EXTRACTION %s\n','NON');
            end
        elseif numel(strfind(line1,'NATURE DES COORDONNEES')) > 0
            if isgeo == 1
                fprintf(fid2,'NATURE DES COORDONNEES\t%s\n','GEOGRAPHIQUES');
            else
                fprintf(fid2,'NATURE DES COORDONNEES\t%s\n','CARTOGRAPHIQUES');
            end
        elseif numel(strfind(line1,'NOMBRE DE LIGNES')) > 0
            fprintf(fid2,'NOMBRE DE LIGNES\t%d\n',nl);
        elseif numel(strfind(line1,'NOMBRE DE COLONNES')) > 0
            fprintf(fid2,'NOMBRE DE COLONNES\t%d\n',nc);
        elseif numel(strfind(line1,'FICHIER BINAIRE')) > 0
            fprintf(fid2,'FICHIER BINAIRE\t%s\n',fi2);
        elseif isext==1 && numel(strfind(line1,'PREMIERE LIGNE EXTRAITE')) > 0
            fprintf(fid2,'PREMIERE LIGNE EXTRAITE\t%d\n',l1);
        elseif isext==1 && numel(strfind(line1,'PREMIERE COLONNE EXTRAITE')) > 0
            fprintf(fid2,'PREMIERE COLONNE EXTRAITE\t%d\n',c1);
        elseif isext==1 && numel(strfind(line1,'LIGNES EXTRAITES')) > 0
            fprintf(fid2,'LIGNES EXTRAITES\t%d\n',ml);
        elseif isext==1 && numel(strfind(line1,'COLONNES EXTRAITES')) > 0
            fprintf(fid2,'COLONNES EXTRAITES\t%d\n',mc);
            % get line numbers, etc for geographic coordinates
        elseif (isgeo==1 && numel(strfind(line1,'LATITUDE DU POINT 1')) > 0)
            fprintf(fid2,'LATITUDE DU POINT 1\t%.10f\n',y1);
        elseif (isgeo==1 && numel(strfind(line1,'LONGITUDE DU POINT 1')) > 0)
            fprintf(fid2,'LONGITUDE DU POINT 1\t%.10ff\n',x1);
        elseif (isgeo==1 && numel(strfind(line1,'PAS LONGITUDE')) > 0)
            fprintf(fid2,'PAS LONGITUDE\t%.10f\n',dc);
        elseif (isgeo==1 && numel(strfind(line1,'PAS LATITUDE')) > 0)
            fprintf(fid2,'PAS LATITUDE\t%.10f\n',dl);
            %             % get line numbers, etc for Cartographic coordinates
            %         elseif (isgeo>1 & numel(strfind(line1,'COORDONNEE Y POINT 0')) > 0)
            %             dums = strread(line1,'%s');
            %             y1 = sscanf(dums{5},'%f');
            %         elseif (isgeo>1 & numel(strfind(line1,'COORDONNEE X POINT 0')) > 0)
            %             dums = strread(line1,'%s');
            %             x1 = sscanf(dums{5},'%f');
            %         elseif (isgeo>1 & numel(strfind(line1,'PAS X')) > 0)
            %             dums = strread(line1,'%s');
            %             dc = sscanf(dums{3},'%f');
            %         elseif (isgeo>1 & numel(strfind(line1,'PAS Y')) > 0)
            %             dums = strread(line1,'%s');
            %             dl = sscanf(dums{3},'%f');
            %         elseif (isgeo>1 & numel(strfind(line1,'LATITUDE ORIGINE EN DEGRES')) > 0)
            %             dums = strread(line1,'%s');
            %             if numel(dums) > 4
            %                 %lat0 = sscanf(dums{5},'%f'); % invokes error
            %                 lat0 = str2double(dums{5}); % returns NAN on error
            %             else
            %                 lat0=NaN;
            %             end
            %         elseif (isgeo>1 & numel(strfind(line1,'LONGITUDE ORIGINE EN DEGRES')) > 0)
            %             dums = strread(line1,'%s');
            %             if numel(dums) > 4
            %                 %lon0 = sscanf(dums{5},'%f');
            %                 lon0 = str2double(dums{5}); % returns NAN on error
            %             else
            %                 lon0=NaN;
            %             end
            %         elseif (isgeo>1 & numel(strfind(line1,'COORDONNEES EN X DE L''ORIGINE')) > 0)
            %             dums = strread(line1,'%s');
            %             if numel(dums) > 5
            %                 %x0 = sscanf(dums{6},'%f');
            %                 x0 = str2double(dums{6}); % returns NAN on error
            %             end
            %         elseif (isgeo>1 & numel(strfind(line1,'COORDONNEES EN Y DE L''ORIGINE')) > 0)
            %             dums = strread(line1,'%s');
            %             if numel(dums) > 5
            %                 %y0 = sscanf(dums{6},'%f')
            %                 y0 = str2double(dums{6}); % returns NAN on error
            %             end
            %         elseif (isgeo>1 & numel(strfind(line1,'HEMISPHERE')) > 0)
            %             dums = strread(line1,'%s');
            %             if numel(dums) > 1
            %                 if strcmpi(dums{2},'NORD')==1
            %                     hemisphere = 'N';
            %                 elseif strcmpi(dums{2},'SUD')==1
            %                     hemisphere = 'S';
            %                 end
            %             else
            %                 if (lat0 > 0)
            %                     hemisphere = 'N';
            %                 elseif (lat0 > 0)
            %                     hemisphere = 'S';
            %                 else
            %                     hemisphere = '?';
            %                 end
            %             end
            %         elseif(isgeo == 2 & numel(strfind(line1,'NUMERO DU FUSEAU')) > 0)
            %             dums = strread(line1,'%s');
            %             if numel(dums) > 3
            %                 %y0 = sscanf(dums{6},'%f')
            %                 iutmzone = str2double(dums{4}); % returns NAN on error
            %             end
        else
            % copy line from input to output
            fprintf(fid2,'%s\n',line1);
        end
    end
end % loop over lines of file

%fi2
%nl
%nc
% if (i == 0)
%     fprintf(1,'ERROR No lines in %s matched key %s\n',igram_list_file,key);
%     error 'No records matched key'
%     return
% end
% fprintf(1,'Read %4d pairs containing key %s\n',i,key);
%
return

