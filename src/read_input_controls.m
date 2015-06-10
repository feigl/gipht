function [xcenter, ycenter, halfwidth, halfheight, npix, pselect, tquake...
   , unitv, ithresh, maxcmd, pixinpatch, maxpix, ianneal, nprocessors, interpcell ...
   , ilist, txtinname, txtoutname, objfun, fitfun, demdescfile, orbfile, cohfile...
   , mpercy, datafilename, nsaruns, parmfilename, saopt6, figopt, printfun, orbopt ...
   , pha2qlsname, phaseprefix, surrogate, verbose]...
   = read_input_controls(fname,runname)
% Read a subregion descriptor file to get the region of interest and pixel selection

% Modifications
% 2009-MAR-25: Lee, version 1
% 2009-MAR-28: Kurt, add loncenter, latcenter, tquake
% 2009-MAR-29: Kurt: add ianneal
% 2009-APR-02: Lee, revised parameter parsing to improve robustness
% 2009-MAY-10: Kurt: allow pselect == 4
% 2009-JUN-16: Kurt: allow pselect == 5
% 2009-JUL-14: Kurt: change name, add txtinname, txtoutname, objfun, fitfun, demdescfile,
% 2009-DEC-08: Kurt: allow pselect = 7
% 2010-JUN-27: Kurt: add maxcmd
% 2010-OCT-05: Kurt: add datafilename
% 2011-APR-02: Kurt: add nsaruns
% 2011-APR-18: Kurt: add parmfilename
% 2011-JUN-14: Kurt: saopt6
% 2011-JUL-03: Kurt: figopt
% 2011-JUL-23: Kurt: maxpix
% 2014-MAR-16: Kurt: surrogate

fprintf(1,'%s begin reading input control file named %s...\n',mfilename,fname);

fid = fopen(fname,'r');
if fid <= 0 
    error(sprintf ('ERROR: Cannot open input control file called %s\n',fname));
end

% read the entire file into idat
%idat=textscan(fid,'%s%f','CommentStyle','%');
%idat=textscan(fid,'%s%s','CommentStyle','%');
%idat=textscan(fid,'%s%s','CommentStyle','%','Delimiter',' ','MultipleDelimsAsOne',1);
%idat=textscan(fid,'%s%s','CommentStyle','%','Delimiter','=','MultipleDelimsAsOne',0,'CollectOutput',1)

% whos idat
% size(idat)
% % 
% keys = idat{:,1}
% strs = idat{:,2}
% 
% 
% nk = numel(keys);	% Number of key words
% nv = numel(strs);	% Number of strings
% % 
% if nk ~= nv
%     error(sprintf ('ERROR: Differing count %d keys and %d values from %s\n',nk,nv,fname));
% end
% 
% for i=1:nk
%    vals{i} = cellstr(sprintf('%s %s',char(keys{i}),char(strs{i})));
% end

k=0;
while 1
   tline = fgetl(fid);
   if ~ischar(tline)
      break
   else
      %i1=findstr(tline,'=');
      i1=strfind(tline,'=');
      if numel(i1) > 1
          error(sprintf('Found more than one equal sign in gipht.in file.\n %s\n',tline));
      end
      %i2=findstr(tline,'%');
      i2=strfind(tline,'%');
      if numel(i2) > 0
         i2 = i2(1);
      else
         i2 = numel(tline);
      end
      if i1 > 0 & i1 < i2
         k=k+1;
         keys{k}=strtrim(tline(1:i1-1));
         vals{k}=tline(i1+1:end);
      end        
   end
end
nv=k;
fclose(fid);
for i=1:nv
   fprintf('%5d %s %s\n',i,keys{i},vals{i});
end
% whos vals
% numel(vals)
% vals


% 
% vals=strcat(keys,vals1)

% for i=1:nv
%    sval=vals1{i};
%    i2=findstr(sval,'%')-1;
%    if i2 == 0
%       i2=numel(sval);
%    end
%    vals{i}=sprintf('%s',sval(1:i2));
%    %vals{i}=textscan(sval,'%s%s','Delimiter','%');
% end
% strmatch('xcenter',keys)
% whos vals{1}
% fprintf(1,'%s\n',char(vals{1}));
% fprintf(1,'%f\n',getval2(keys,vals,'xcenter','f',nv);


% Load the values associated with keys
%  If a key isn't found it will be defaulted later
xcenter      = getval2(keys,vals,'xcenter','f',nv);
ycenter      = getval2(keys,vals,'ycenter','f',nv);
halfwidth    = getval2(keys,vals,'halfwidth','f',nv);
halfheight   = getval2(keys,vals,'halfheight','f',nv);
npix         = getval2(keys,vals,'npix','f',nv);
pselect      = getval2(keys,vals,'pselect','f',nv);
unitv_east   = getval2(keys,vals,'unitv_east','f',nv);
unitv_north  = getval2(keys,vals,'unitv_north','f',nv);
unitv_up     = getval2(keys,vals,'unitv_up','f',nv);
tquake       = getval2(keys,vals,'tquake','f',nv);
ithresh      = getval2(keys,vals,'ithresh','f',nv);
ianneal      = getval2(keys,vals,'anneal','f',nv);
pixinpatch   = getval2(keys,vals,'pixinpatch','f',nv);
maxpix       = getval2(keys,vals,'maxpix','f',nv);
nprocessors  = getval2(keys,vals,'nprocessors','f',nv);
interpcell   = getval2(keys,vals,'interpcell','f',nv);
mpercy       = getval2(keys,vals,'mpercy','f',nv);
maxcmd       = getval2(keys,vals,'maxcmd','f',nv);
nsaruns      = getval2(keys,vals,'nsaruns','f',nv);
saopt6       = getval2(keys,vals,'saopt6','f',nv);
%figopt       = getval2(keys,vals,'figopt','f',nv);
figopt       = getval2(keys,vals,'figopt','b',nv);
orbopt       = getval2(keys,vals,'orbopt','f',nv);
surrogate    = getval2(keys,vals,'surrogate','f',nv);
verbose      = getval2(keys,vals,'verbose','f',nv);

% New string items 2009-JUL-14
objfun       = getval2(keys,vals,'objfun','s',nv);
fitfun       = getval2(keys,vals,'fitfun','s',nv);
txtinname    = getval2(keys,vals,'txtinname','s',nv);
txtoutname   = getval2(keys,vals,'txtoutname','s',nv);
demdescfile  = getval2(keys,vals,'demdescfile','s',nv);
ilist        = getval2(keys,vals,'ilist','s',nv);
orbfile      = getval2(keys,vals,'orbfile','s',nv);
cohfile      = getval2(keys,vals,'cohfile','s',nv);
datafilename = getval2(keys,vals,'datafilename','s',nv);
parmfilename = getval2(keys,vals,'parmfilename','s',nv);
printfun     = getval2(keys,vals,'printfun','s',nv);
pha2qlsname  = getval2(keys,vals,'pha2qlsname','s',nv);
phaseprefix  = getval2(keys,vals,'phaseprefix','s',nv);


% Make sure a value is found
%  this code can catch typos in keys
%  or assign default values for omitted entries.
if numel(xcenter) == 0
    warning(sprintf ('WARNING: xcenter not yet defined in %s\n',fname));
end
if numel(ycenter) == 0
    warning(sprintf ('WARNING: ycenter not yet defined in %s\n',fname));
end
if numel(halfwidth) == 0
	halfwidth = 120;	% default
end
if numel(halfheight) == 0
	halfheight = 100;	% default
end
if numel(npix) == 0
	npix = 199;	% default
end
if numel(pselect) == 0
	pselect = 1; % default to random
end
if numel(unitv_east) == 0
	unitv_east = 0.4;	% default
end
if numel(unitv_north) == 0
	unitv_north = -0.1;	% default
end
if numel(unitv_up) == 0
	unitv_up = 0.9;	% default
end
if numel(tquake) == 0
	tquake = 0;	% default
end
if numel(ithresh) == 0
   ithresh = 8; % threshold circular mean deviation misfit from mean direction in DN [0 127]
end
if numel(maxcmd) == 0
   maxcmd = 16; % threshold circular mean deviation from ramp in DN [0 127]
end
if numel(pixinpatch) == 0
   pixinpatch = 16;   % mininum number of good pixels in a quadtree patch
end
if numel(maxpix) == 0
   maxpix = 16;   % maximum number of pixels in a quadtree patch
end
if numel(ianneal) == 0
   %ianneal = 1;	% default
   ianneal = 0;	% new default 20130113
end
if numel(objfun) == 0
	objfun = 'funcostrarcm'; % mininum angle, assumes zero mean, using matlab function;	% default
end
if numel(fitfun) == 0
   fitfun ='funseparable20';   % 
end
if numel(printfun) == 0
   printfun = 'printnull';   % 
end
if numel(txtinname) == 0
   txtinname ='parameters.gin'; % new format
end
if numel(txtoutname) == 0
   txtoutname=sprintf('%s_%s',runname,strrep(txtinname,'.gin','.out'));
end
if numel(demdescfile) == 0
    demdescfile = 'dem_descriptor.dat'; 
end
if numel(ilist) == 0
   ilist = 'file_names.dat';
end
if numel(nprocessors) == 0
   nprocessors = 1;
end
if numel(interpcell) == 0
   interpcell = 1000; % dimension of cell (in same units as xcenter and ycenter [meters or degrees]) for interpolation in gipht_step4
end
if numel(mpercy) == 0
   mpercy = 0.0284; % assume ERS or ENVISAT
end
if numel(nsaruns) == 0
   nsaruns = 1; % 1 run for SA
end
if numel(saopt6) == 0
   saopt6 = 0; % do not follow anneal4 with optimization
end
if numel(figopt) == 0
   figopt = 0; % no fancy plots
end
if numel(orbopt) == 0
   orbopt = 0; % do not handle estimate orbital parameters
end
if numel(pha2qlsname) == 0
    pha2qlsname ='pha2qls';   % 
end
if numel(phaseprefix) == 0
    phaseprefix ='psp';   % 
end
if numel(surrogate) == 0
   surrogate = 0; % use exact (rather than approximate) version of fitting function
end
if numel(verbose) == 0
   verbose = 0; % be quiet
end

% orbit file
if numel(orbfile) == 0
   %orbfile = 'master.orb';
   fprintf(1,'WARNING: No orbit file named in %s. Assuming constant look vector.\n',fname);
   % check unit vector
   unitv = [unitv_east unitv_north unitv_up];
   for i=1:3
      if abs(unitv(i)) > 1.0
         error(sprintf('ERROR: %ith component of UNITV has invalid value %f in %s\n'...
            ,i,unitv(i),fname));
      end
   end
   i=3;
%    if unitv(i) < 0
%       error(sprintf('ERROR: %ith component of UNITV has invalid value %f in %s\n'...
%          ,i,unitv(i),fname));
%    end
   unorm = abs(norm(unitv) - 1.0);
   if unorm > 1e-2
      error(sprintf('ERROR: norm of UNITV differs from unity by %f in %s\n',unorm,fname));
   end
else
   unitv = [NaN NaN NaN];
   if orbfile(end) == '/'
      fprintf(1,'Will look for individual orbit files in directory named %s\n',orbfile)
   else
      if fexist(orbfile) <= 0
         error(sprintf('Cannot find file %s as listed in %s\n',orbfile,fname));
      end
   end
end
% coherence file
if pselect == 3
   if numel(cohfile) == 0
      fprintf(1,'WARNING: No coherence file named in %s\n',fname);
      %cohfile = 'coh.byt';
   else
      if fexist(cohfile) <= 0
         error(sprintf('Cannot find file %s as listed in %s\n',cohfile,fname));
      end
   end
end


% Check that files exist
Cfiles = {txtinname demdescfile ilist parmfilename};
for i=1:numel(Cfiles);
    if numel(Cfiles{i}) > 0
        fprintf('%s\n',Cfiles{i});
        if fexist(Cfiles{i}) <= 0
            error(sprintf('Cannot find file %s as listed in %s\n',Cfiles{i},fname));
        end
    end
end
% Check that files exist
Cfiles = {datafilename};
for i=1:numel(Cfiles);
    if numel(Cfiles{i}) > 0
        fprintf('%s\n',Cfiles{i});
        if fexist(Cfiles{i}) <= 0
            warning(sprintf('Cannot find file %s as listed in %s\n',Cfiles{i},fname));
        end
    end
end

% check that functions exist
% if exist('which') == 5
%     NOTE: COMMAND "WHICH" DOES NOT APPEAR TO WORK IN COMPILED VERSION
%     cmds = {fitfun,printfun,objfun}
%     for i=1:numel(cmds)
%         cmd1 = char(cmds{i});
%         fprintf(1,'Using "which" command to search for function named %s\n',cmd1);
%         result = which(cmd1)
%         if numel(result) > 0
%             fprintf(1,'Found fitting function named %s %s\n',cmd1,result);
%         else
%             error(sprintf('Command "which" failed while looking for function named %s\n',cmd1));
%         end
%     end
% else
%     warning('Cannot find which command');
% end

% check that functions exist
cmds = {fitfun,printfun,objfun}
for i=1:numel(cmds)
    cmd1 = char(cmds{i});
    fprintf(1,'Using "exist" command to search for function named %s\n',cmd1);
    result = exist(cmd1,'file');
    if result > 0
        fprintf(1,'Found function named %s %d\n',cmd1,result);
    else
        error(sprintf('Command "which" failed while looking for function named %s\n',cmd1));
    end
end

if numel(parmfilename) == 0
   parmfilename = 'gipht_parm.pst';
end

% Perform sanity checks
if  ~ismember(pselect,[0,1,2,3,5,7,9]);
    error(sprintf ('ERROR: PSELECT invalid value in %s\n',fname));
end

if  ianneal < 0  || ianneal > 6 
    error(sprintf ('ERROR: ANNEAL invalid value in %s\n',fname));
end

if  nsaruns <= 0  || nsaruns > 10 
    error(sprintf ('ERROR: NSARUN invalid value in %s\n',fname));
end

if ithresh < 0 || ithresh > 127
    error(sprintf('ERROR: ithresh has invalid value of %f in %s\n',ithresh,fname));
end

if maxcmd < 0 || maxcmd > 127
    error(sprintf('ERROR: maxcmd has invalid value of %f in %s\n',maxcmd,fname));
end

if pixinpatch < 0 || pixinpatch > 100
    error(sprintf('ERROR: pixinpatch has invalid value of %f in %s\n',pixinpatch,fname));
end

if interpcell < 0 || interpcell > 100e3
    error(sprintf('ERROR: interpcell has invalid value of %f in %s\n',interpcell,fname));
end

if pselect < 3
    if npix < 10 || npix > (2*halfwidth + 1) * (2*halfheight +1)
        error(sprintf ('ERROR: NPIX is < 10 or > H*W in %s\n',fname));
    end
end

if figopt < 0 || figopt > 8
    error(sprintf ('ERROR: FIGOPT invalid value %d in %s\n',figopt,fname));
end
% % vector displacement requires interpolation
% if bitget(figopt,3) == 1 && bitget(figopt,2) == 0
%     warning(sprintf('vector displacement requires interpolation\n'));
%     figopt = bitset(figopt,2,1);
%     fprintf(1,'figopt = %4s\n',dec2bin(figopt,4));
% end

if  ~ismember(orbopt,[0,1]);
    error(sprintf('ERROR: orbopt has invalid value of %f in %s\n',orbopt,fname));
end

% if nprocessors > 0
%     if numel(which('matlabpool')) > 0
%         matlabpoolsize = matlabpool('size');
%         if matlabpoolsize > 0
%            try 
%                matlabpool close force
%            catch
%                warning('issues with matlab pool');
%            end
%         else
%             warning('matlabpoolsize is %d\n',matlabpoolsize);
%         end
%     else
%         error(sprintf('Request is for nprocessors = %d BUT matlab distributed tool kit not installed.\n',nprocessors));
%     end
% end
% 
% In Matlab version 8.5.0.197613 (R2015a)
%  matlabpool has been removed. Use PARPOOL instead.

if nprocessors > 1 
    if numel(which('parpool')) > 0
        current_parpool = gcp('nocreate');
        if numel(current_parpool) > 0
            warning('Deleting current pool...');
            delete(gcp);
        end
        %     pool = parpool(numWorkers) creates and returns a pool with the
        %     specified number of workers. numWorkers can be a positive integer or
        %     a range specified as a 2-element vector of integers. If numWorkers is
        %     a range, the resulting pool has size as large as possible in the
        %     range requested.
        current_parpool= parpool([1,nprocessors]);
    else
        error(sprintf('Request is for nprocessors = %d BUT matlab distributed tool kit not installed.\n',nprocessors));
    end
end

if mpercy < 0.001 || mpercy > 1.000
    error(sprintf('ERROR: mpercy has invalid value of %f in %s\n',mpercy,fname));
end

allowed = {'pha','psp','cln','psm'};
if ismember(phaseprefix,allowed) ~= 1
    allowed
    error(sprintf('ERROR: phaseprefix %s is not on list of allowed values above.\n',phaseprefix));
end

fprintf(1,'%s completed reading input controls.\n',mfilename);
return
