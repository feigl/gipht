function grd_file_names = aria_to_gmt_grd(dirname,input_file_name)
%% convert files from ARIA to GMT grd format for use with GIPhT
% 20180926 Kurt Feigl

%% parse its xml file
 
XML = xml2struct(strcat(dirname,filesep,input_file_name,'.xml'));

% XML.imageFile
% XML.imageFile.property
% XML.imageFile.component
% XML.imageFile.component(1)
% XML.imageFile.component{1}
% XML.imageFile.component{1}.property
% XML.imageFile.component{1}.property.{1}
% XML.imageFile.component{1}.property{1}
% XML.imageFile.component{1}.Attributes
% XML.imageFile.component{:}.Attributes

% parse properties
for i=1:numel(XML.imageFile.property)
    fileProperty1 = XML.imageFile.property{i};
    textvalue1  = fileProperty1.value.Text;
    name1 = fileProperty1.Attributes.name;
    fprintf(1,'%20s : %20s\n',name1,textvalue1);
    dvalue1 = str2double(textvalue1);
    if isfinite(dvalue1) == 1
        META_PAIR.(name1)=dvalue1;
    else
        META_PAIR.(name1)=textvalue1;
    end
end
%META_PAIR


%% now parse components for pair
components = XML.imageFile.component;
ncomponents = numel(XML.imageFile.component);
for i=1:ncomponents
    nproperties = numel(XML.imageFile.component{i}.property);
    doc1 = components{i}.doc.Text;
    if strcmp(components{i}.factorymodule.Text, 'isceobj.Image') == 1
        for j=1:nproperties
            name1      = XML.imageFile.component{i}.property{j}.Attributes.name;
            textvalue1 = XML.imageFile.component{i}.property{j}.value.Text;
            fprintf(1,'%45s %20s : %20s\n',doc1,name1,textvalue1);
            dvalue1 = str2double(textvalue1);
            if isfinite(dvalue1) == 1
                switch i
                    case 1
                        COORD1.(name1)=dvalue1;
                        COORD1.doc    =doc1;
                    case 2
                        COORD2.(name1)=dvalue1;
                        COORD2.doc    =doc1;
                    otherwise
                        error(sprintf('unknown compnonent %d',i));
                end
            end
        end
    end
end
% COORD1
% COORD2

%% TODO: parse VRT file instead
% From: Kurt Feigl <feigl@wisc.edu> Date: Friday, September 28, 2018 at
% 3:31 PM To: "Bekaert, David (334H)" <David.Bekaert@jpl.nasa.gov> Subject:
% Re: [aria-help] Re: Description of ARIA interferogram output files?
%  
% Hi David,
%  
% Thanks for your rapid reply. It would be great if these descriptions were
% included in the XML file.
%  
% Along the same lines, I found an error in the attached file. The ?ending
% value? for the coordinate is not a longitude.
%  
% Kurt
%  
% On  2018-Sep-29, at 12:45 , Bekaert, David (334H) <David.Bekaert@jpl.nasa.gov> wrote:
% Looks like something erroneous. Better to use the vrt information, which
% is accurate.
%  
% The xml information is specific to isce and currently only kept to be
% backward compatible. I suspect this leaked through as the ending
% information never gets used (i.e. directly laoded from the xml) Typically
% get calculated on the fly from starting value, samples and delta to get
% to the ending value I forwarded this to the ISCE team.
%  
% Cheers, D.
%  

dlat = COORD2.delta
dlon = COORD1.delta

nlat = COORD2.size
nlon = COORD1.size

lat1 = COORD2.startingvalue;
lon1 = COORD1.startingvalue;

lat2 = lat1 + nlat*dlat;
lon2 = lon1 + nlon*dlon;


%% check sizes
if META_PAIR.width == COORD1.size
    mcols = META_PAIR.width
else
    error('miscount in width')
end
if META_PAIR.length == COORD2.size
    nrows = META_PAIR.length
else
    error('miscount in length')
end

%% check size of file
DIRDATA = dir(strcat(dirname,filesep,input_file_name));
nbytes = DIRDATA.bytes

% if META_PAIR.width * META_PAIR.length * 8 ~= nbytes
%     error('miscount in number of bytes');
% end
npix = META_PAIR.width * META_PAIR.length
if npix*8 == nbytes
    nfiles = 2
    nbands = 2
elseif npix*4 == nbytes
    nfiles = 1
    nbands = 1
else
    error('miscount in number of bytes');
end


switch input_file_name
    case {'topophase.flat.geo','filt_topophase.flat.geo'}
        %% read the phase file
        z = read_cr4(strcat(dirname,filesep,input_file_name),META_PAIR.length,META_PAIR.width);
        %nfiles = 2;       
        field{1} = angle(z);
        field{2} = abs(z);
        % mask phase where amplitude is zero
        %ibad = find(abs(field{2}) < eps);
        % less than 5%
        ibad = find(abs(field{2}) < quantile(colvec(abs(field{2})),0.05));
        field{1}(ibad) = nan;
        clear z;
        description{1} =  'Phase     [radians]';tag{1} = 'pha';
        description{2} =  'Amplitude [dimless]';tag{2} = 'amp';
        
    case 'los.rdr.geo'        
        %% Get pointing vector from target on ground to radar sensor on sattelite
        % look in file named merged/los.rdr.geo.xml
        % %    <property name="description">
        %         <value>["['Two channel Line-Of-Sight geometry image (all angles
        %         in degrees). Represents vector drawn from target to platform. \\n
        %         Channel 1: Incidence angle measured from vertical at target (always +ve).\\n
        %         Channel 2: Azimuth angle measured from North in Anti-clockwise direction.']"]</value>
        %         <doc>Image description</doc>
        %nbands = 2;
        r = multibandread(strcat(dirname,filesep,input_file_name), [nrows mcols nbands], ...
            'float32', 0, 'bil', 'ieee-le');
        %nfiles = 2;       
        field{1} = r(:,:,1);
        field{2} = r(:,:,2);
        clear r;
        description{1} =  'Incidence angle measured from vertical at target (always +ve). [degrees]';tag{1}='inc';
        description{2} =  'Azimuth angle measured from North in Anti-clockwise direction  [degrees]';tag{2}='azi';            
    case 'phsig.cor.geo'
        %% read PHAse SIGma?
        
        %
        %
        % On  2018-Sep-28, at 16:22 , Bekaert, David (334H)
        % <David.Bekaert@jpl.nasa.gov> wrote:
        %
        % Hi,
        %
        % phsig.cor.geo is a single band file containing the correlation (value
        % between 0-1 and dimensionless) of the filtered interferogram
        % filt_topophase.unw.geo which you can use as proxi for coherence. If you
        % decide to use the unfiltered unwrapped interferogram then you should use
        % the .cor file which is a two band file with the coherence in the second
        % band, again between 0-1 and dimensionless. The filtered interferogram
        % filt_topophase.unw.geo is a two band file with the first band the
        % amplitude of the interferogram and the second band the unwrapped and
        % filtered phase in radians.
        %
        % D.
        %file_name_cor = strcat(dirname,filesep,'merged',filesep,'phsig.cor.geo')
        
        %nbands = 1;
        r = multibandread(strcat(dirname,filesep,input_file_name), [nrows mcols nbands], ...
            'float32', 0, 'bil', 'ieee-le');
        %SG1 = r(:,:,1);
        nfiles = 1;       
        field{1} = r(:,:,1);
        clear r;
        description{1} =  'phsig.cor.geo what is this? [??]';tag{1} = 'sig';        
    case 'filt_topophase.unw.geo'              
        %% unwrapped
        %nbands = 2;
        r = multibandread(strcat(dirname,filesep,input_file_name), [nrows mcols nbands], ...
            'float32', 0, 'bil', 'ieee-le');
        %nfiles = 2;       
        field{1} = r(:,:,1);
        field{2} = r(:,:,2);
        clear r;
        description{1} =  'Amplitude        [dimless]';tag{1} = 'uw1';
        description{2} =  'Unwrapped phase  [radians]';tag{2} = 'uw2';                   
    otherwise
        error(sprintf('unknown type of file %s\n',input_file_name));
end

%% flip if necessary
if  dlat < 0
    warning('Flipping in latitude');
    field{1} = flipud(field{1});
    if nfiles == 2
        field{2} = flipud(field{2});
    end
    latax = lat2:abs(dlat):lat1;
else
    latax = lat1:abs(dlat):lat2;
end

%% flop if necessary
if dlon < 0
    warning('flopping in longitude');
    field{1} = fliplr(field{1});
    if nfiles == 2
        field{2} = fliplr(field{2});
    end
    lonax = lon2:abs(dlon):lon1;
else
    lonax = lon1:dlon:lon2;
end

% build a header
for i=1:nfiles
    %     grdfilename = sprintf('%s_%3s.grd',input_file_name,char(field_name1))
    %     INFO(i).title =  grdfilename;
    INFO(i).conventions =  'COARDS/CF-1.0';
    INFO(i).gmtversion =  '4.x';
    INFO(i).command =  sprintf('File written by %s %\n',mfilename,datestr(now));
    INFO(i).ispixelreg =  1;
    INFO(i).dx =  abs(dlon);
    INFO(i).dy =  abs(dlat);
    INFO(i).xmin =  nanmin(lonax);
    INFO(i).xmax =  nanmax(lonax);
    INFO(i).ymin =  nanmin(latax);
    INFO(i).ymax =  nanmax(latax);
    switch i
        case 1
            INFO(i).zmin =  nanmin(colvec(field{1}));
            INFO(i).zmax =  nanmax(colvec(field{1}));
        case 2
            INFO(i).zmin =  nanmin(colvec(field{2}));
            INFO(i).zmax =  nanmax(colvec(field{2}));
        otherwise
            error(sprintf('bad i %d',i));
    end
    
    INFO(i).nx =    mcols;
    INFO(i).ny =    nrows;
    INFO(i).xname =  'Longitude in degrees';
    INFO(i).yname =  'Latitude  in degrees';
    INFO(i).zname =  description{i};
end


%% Write a GMT grid file
for i=1:nfiles
    grd_file_names{i} = strcat(dirname,filesep,input_file_name,'.',tag{i},'.grd');
    grdwrite3(lonax,latax,field{i},grd_file_names{i},INFO(i));
    fprintf(1,'Wrote %s\n',grd_file_names{i});
    dirlist = dir(grd_file_names{i})
    
%     %% Make a map of it
%     if contains(grd_file_names{i},'AMP') == 1
%         cmap = gray;
%     else
%         cmap = jet;
%     end
%     H=map_grd(grd_file_names{i},cmap);
end


return

end

