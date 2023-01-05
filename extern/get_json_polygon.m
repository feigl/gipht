function [X,Y] = get_json_polygon(fileName)

%https://www.mathworks.com/matlabcentral/answers/474980-extract-info-from-json-file-by-matlab
%fileName = 'filename.json'; % filename in JSON extension
fid = fopen(fileName); % Opening the file
raw = fread(fid,inf); % Reading the contents
str = char(raw'); % Transformation
fclose(fid); % Closing the file
XY = jsondecode(str); 
X = XY.features.geometry.coordinates(:,:,1);
Y = XY.features.geometry.coordinates(:,:,2);
end