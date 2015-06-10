function npatches = psfilt3(pha_in,amp_in,pha_out,ncol,alpha)
%function npatches = psfilt3(pha_in,amp_in,alpha,ncols,pha_out)
% Apply power spectral filtering using Goldstein and Werner algorithm
% pha_in  == name of input phase file in Diapason .pha format
% amp_in  == name of input amplitude file in Diapason .byt format
% alpha   == exponent
% ncols   == number of columns
% pha_out == name of output phase file in Diapason .pha format
% 20140531 Kurt Feigl

if nargin == 5
    
    npatches = NaN;
    
    cmd1 = get_executable_name('ps_filt3.c');
%   usage: ./ps_filt3.maci64 <pha.oct> <amp.oct> <smp.oct> <ncol> <alpha> 
% 
% input parameters: 
% pha.oct  (Input)  Phase     (Unsigned 8Bit Integer -128 = -pi; 127 = pi)
% smp.oct  (Output) Smoothed Phase (Unsigned 8Bit Integer -128 = -pi; 127 = pi)
% ncols    (Input)  Number of columns in all 3 files (pixels per line);
% alpha    (Input)  Alpha exponent in Goldstein & Werner technique (recommend: 0.0 < 0.7 < 1.0);
  
    % file names and sizes, ithresh = 4 DN, minpix = 16
    cmd2=sprintf('%s %s %s %d %f'...
        ,pha_in,amp_in,pha_out,ncol,alpha)
    % complete command line
    cmd3=sprintf('%s %s',cmd1,cmd2);
    fprintf(1,'Starting pha2qls with command line:\n%s\n',cmd3);
    %[unixstat,unixout] = unix(cmd3);
    [unixstat,unixout] = system(cmd3);
    if unixstat == 0
        fprintf(1,'%s successful.\n',cmd1);
        %whos unixout
        %unixout
        key = 'N(OK patches)  =';
        k = strfind(unixout,key)+numel(key)+1;
        if k > 0
            npatches = str2num(unixout(k:k+6));
        end
    else
        error(['Error executing psfilt3 ' unixout]);
    end
else  
    error('Incorrect number of arguments\n');
end
return
    
    
    
    
