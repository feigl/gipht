function attrvalue = getAtt(ncid,varid,attname,output_datatype)
%netcdf.getAtt Return netCDF attribute.
%   attrvalue = netcdf.getAtt(ncid,varid,attname) reads an attribute
%   value.  The class of attrvalue will match that of the internal
%   attribute datatype.  For example, if the attribute has netCDF 
%   datatype NC_INT, then the class of the output data will be int32.
%   If an attribute has netCDF datatype NC_BYTE, it will result in an
%   int8 value.
%
%   attrValue = netcdf.getAtt(ncid,varid,attname,attClass) specifies the
%   class of attrValue in the final input argument.  This cannot be used to
%   convert between character and numeric attributes.
%
%   The list of allowable datatype strings consists of 'double', 
%   'single', 'uint64', 'int64', 'uint32', 'int32', 'uint16', 'int16', 
%   'int8', and 'uint8'.
%   
%   Example:
%       ncid = netcdf.open('example.nc','NOWRITE');
%       varid = netcdf.inqVarID(ncid,'temperature');
%       add_offset = netcdf.getAtt(ncid,varid,'add_offset','single');
%       netcdf.close(ncid);
%
%   This function corresponds to the "nc_get_att" family of functions in 
%   the netCDF library C API.
%
%   Please read the files netcdfcopyright.txt and mexnccopyright.txt for 
%   more information.
%
%   See also netcdf, netcdf.putAtt.
%

%   Copyright 2008-2013 The MathWorks, Inc.

if nargin > 2
    attname = convertStringsToChars(attname);
end

if nargin > 3
    output_datatype = convertStringsToChars(output_datatype);
end

persistent nc_classes;
nc_classes = { 'double','float','single','int64','uint64', ...
               'int','int32','uint','uint32', 'short','int16', ...
               'ushort','uint16', 'schar','int8','uchar', ...
               'uint8','char','text'};
switch ( nargin )
    case 3
        % Use the internal datatype to determine the output class.
        fprintf(1,'attname is %s\n',attname);
        xtype = netcdflib('inqAtt',ncid,varid,attname)
        switch ( xtype )
            case 11
                funcstr = 'getAttUint64';
            case 10
                funcstr = 'getAttInt64';
            case 9
                funcstr = 'getAttUint';
            case 8
                funcstr = 'getAttUshort';
            case 7
                funcstr = 'getAttUbyte';
            case 6
                funcstr = 'getAttDouble';
            case 5
                funcstr = 'getAttFloat';
            case 4
                funcstr = 'getAttInt';
            case 3
                funcstr = 'getAttShort';
            case 2
                funcstr = 'getAttText';
            case 1
                funcstr = 'getAttSchar';
            otherwise
                %error(message('MATLAB:imagesci:netcdf:unhandledAttributeDatatype', xtype));
                warning(message('MATLAB:imagesci:netcdf:unhandledAttributeDatatype', xtype));
                funcstr = [];
        end
        
    case 4
        % In this case we determine the funcstr from the specified class.
        output_datatype = validatestring(output_datatype,nc_classes);
        switch ( output_datatype )
            case { 'double' }
                funcstr = 'getAttDouble';
            case { 'float', 'single' }
                funcstr = 'getAttFloat';
            case { 'int64' }
                funcstr = 'getAttInt64';
            case { 'uint64' }
                funcstr = 'getAttUint64';
            case { 'int', 'int32' }
                funcstr = 'getAttInt';
            case { 'uint', 'uint32' }
                funcstr = 'getAttUint';
            case { 'short', 'int16' }
                funcstr = 'getAttShort';
            case { 'ushort', 'uint16' }
                funcstr = 'getAttUshort';
            case { 'schar', 'int8' }
                funcstr = 'getAttSchar';
            case { 'uchar', 'uint8' }
                funcstr = 'getAttUchar';
            case { 'char', 'text' }
                funcstr = 'getAttText';
        end
        
end
        
        
attrvalue = netcdflib(funcstr, ncid, varid, attname);            
