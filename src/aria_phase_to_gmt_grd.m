function [AMP,PHA,INC,AZI,SG1,SG2,UW1,UW2,COORD1,COORD2,META_MASTER,META_SLAVE,META_PAIR] = read_aria_interferogram(dirname)


%% read xml file for master 
file_name_master = strcat(dirname,filesep,'master');
XML_MASTER = xml2struct(strcat(file_name_master,'.xml'))
% XML_MASTER.productmanager_name
% XML_MASTER.productmanager_name.property
% XML_MASTER.productmanager_name.property(1)
% XML_MASTER.productmanager_name.component
% XML_MASTER.productmanager_name.component.property
% XML_MASTER.productmanager_name.component.Attributes
% XML_MASTER.productmanager_name.component.property{1}
% XML_MASTER.productmanager_name.component.property{1}.value
% XML_MASTER.productmanager_name.component.component.component{1}.component{1}.property{1}.Attributes


% parse properties
% for i=1:numel(XML_M.productmanager_name.component.property)
%     productProperty1 = XML_M.productmanager_name.component.property{i};
for i=1:numel(XML_MASTER.productmanager_name.component.component.component)
    component1 = XML_MASTER.productmanager_name.component.component.component{i};
    name1 = component1.factoryname.Text;
    for j=1:numel(component1.property)
        name2       = component1.property{j}.Attributes.name;
        doc2        = component1.property{j}.doc.Text;      
        textvalue2  = component1.property{j}.value.Text;
        fprintf(1,'%10s : %20s : %40s %40s\n',name1,name2,textvalue2,doc2);
        dvalue2 = str2double(textvalue2);
        if isfinite(dvalue2) == 1
            value = dvalue2;
        else
            value = textvalue2;
        end
        %META_M.(name2)=struct('value',value,'doc',doc2);
        %META_M.(name2) = {value,doc2};
        META_MASTER.(name2).value = value;
        META_MASTER.(name2).doc   = doc2;
    end
end

%% read xml file for slave
file_name_slave = strcat(dirname,filesep,'master');
XML_SLAVE = xml2struct(strcat(file_name_slave,'.xml'))
% XML_SLAVE.productmanager_name
% XML_SLAVE.productmanager_name.property
% XML_SLAVE.productmanager_name.property(1)
% XML_SLAVE.productmanager_name.component
% XML_SLAVE.productmanager_name.component.property
% XML_SLAVE.productmanager_name.component.Attributes
% XML_SLAVE.productmanager_name.component.property{1}
% XML_SLAVE.productmanager_name.component.property{1}.value
% XML_SLAVE.productmanager_name.component.component.component{1}.component{1}.property{1}.Attributes


% parse properties
% for i=1:numel(XML_M.productmanager_name.component.property)
%     productProperty1 = XML_M.productmanager_name.component.property{i};
for i=1:numel(XML_SLAVE.productmanager_name.component.component.component)
    component1 = XML_SLAVE.productmanager_name.component.component.component{i};
    name1 = component1.factoryname.Text;
    for j=1:numel(component1.property)
        name2       = component1.property{j}.Attributes.name;
        doc2        = component1.property{j}.doc.Text;      
        textvalue2  = component1.property{j}.value.Text;
        fprintf(1,'%10s : %20s : %40s %40s\n',name1,name2,textvalue2,doc2);
        dvalue2 = str2double(textvalue2);
        if isfinite(dvalue2) == 1
            value = dvalue2;
        else
            value = textvalue2;
        end
        %META_M.(name2)=struct('value',value,'doc',doc2);
        %META_M.(name2) = {value,doc2};
        META_SLAVE.(name2).value = value;
        META_SLAVE.(name2).doc   = doc2;
    end
end




%% phase file for pair
%file_name_pha = strcat(dirname,filesep,'merged',filesep,'topophase.flat.geo');
file_name_pha = strcat(dirname,filesep,'merged',filesep,'filt_topophase.flat.geo');

% parse its xml file
XML = xml2struct(strcat(file_name_pha,'.xml'));

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
    %     switch(lower(name1))
    %         case 'length'
    %             length = sscanf(char(textvalue1),'%g',1);
    %         case 'width'
    %             width = sscanf(char(textvalue1),'%g',1);
    % %         otherwise
    % %             error(sprintf('unknown attrib1 %s',attrib1));
    %     end
    dvalue1 = str2double(textvalue1);
    if isfinite(dvalue1) == 1
       META_PAIR.(name1)=dvalue1;
    else
       META_PAIR.(name1)=textvalue1;
    end 
end
META_PAIR


%% now parse components for pair
components = XML.imageFile.component;
ncomponents = numel(XML.imageFile.component);
for i=1:ncomponents
    nproperties = numel(XML.imageFile.component{i}.property);
    doc1 = components{i}.doc.Text;
    if strcmp(components{i}.factorymodule.Text, 'isceobj.Image') == 1
        for j=1:nproperties
%             property1  = XML.imageFile.component{i}.property{j};
%             attribute1 = XML.imageFile.component{i}.property{j}.Attributes;
            name1      = XML.imageFile.component{i}.property{j}.Attributes.name;
            textvalue1 = XML.imageFile.component{i}.property{j}.value.Text;
            fprintf(1,'%45s %20s : %20s\n',doc1,name1,textvalue1);
            %             if contains(textvalue1,'.') == true
            %                 dvalue1 = sscanf(char(textvalue1),'%g',1);
            %               if isnumeric(dvalue1) == true
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
DIRDATA = dir(file_name_pha);
nbytes = DIRDATA.bytes

if META_PAIR.width * META_PAIR.length * 8 ~= nbytes
    error('miscount in number of bytes');
end


%% read the phase file
%z = read_cr4(file_name_pha,nrows,mcols);
z = read_cr4(file_name_pha,META_PAIR.length,META_PAIR.width);
PHA = angle(z);
AMP = abs(z);

whos PHA
whos AMP

%% Get pointing vector from target on ground to radar sensor on sattelite
% look in file named merged/los.rdr.geo.xml
% %    <property name="description">
%         <value>["['Two channel Line-Of-Sight geometry image (all angles
%         in degrees). Represents vector drawn from target to platform. \\n
%         Channel 1: Incidence angle measured from vertical at target (always +ve).\\n                
%         Channel 2: Azimuth angle measured from North in Anti-clockwise direction.']"]</value>
%         <doc>Image description</doc>
% read the file
% nbytes1 = 301710136 % from ls
% nbytes2 = 2 * 4 * nrows * mcols
% if nbytes1 ~= nbytes2
%     error('miscount in los file');
% end
file_name_los = strcat(dirname,filesep,'merged',filesep,'los.rdr.geo');
%[INC,AZI] = read_bil2r4(file_name_los,nrows,mcols);
nbands = 2;
r = multibandread(file_name_los, [nrows mcols nbands], ...
                     'float32', 0, 'bil', 'ieee-le'); 
INC = r(:,:,1);
AZI = r(:,:,2);
clear r;

%% read PHAse SIGma?
file_name_cor = strcat(dirname,filesep,'merged',filesep,'phsig.cor.geo')
nbands = 1;
r = multibandread(file_name_cor, [nrows mcols nbands], ...
    'float32', 0, 'bil', 'ieee-le');
SG1 = r(:,:,1);
SG2 = nan(size(SG1));


%% unwrapped
file_name_unw = strcat(dirname,filesep,'merged',filesep,'filt_topophase.unw.geo');
nbands = 2;
r = multibandread(file_name_unw, [nrows mcols nbands], ...
    'float32', 0, 'bil', 'ieee-le');
UW1 = r(:,:,1);
UW2 = r(:,:,2);

whos UW1
whos UW2

return

end

