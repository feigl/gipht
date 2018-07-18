 function OPT = read_input_controls(OPT)
% function OPT = read_input_controls(OPT,OPT.fning)
    % function [xcenter, ycenter, halfwidth, halfheight, npix, pselect, tquake...
%    , unitv, ithresh, maxcmd, pixinpatch, maxpix, ianneal, nprocessors, interpcell ...
%    , ilist, fnparin, fnparout, objfun, fitfun, demdescfile, orbfile, cohfile...
%    , mpercy, datafilename, nsaruns, parmfilename, saopt6, figopt, printfun, orbopt ...
%    , pha2qlsname, phaseprefix, surrogate, verbose]...

% read control parameters from file named OPT.fning (e.g., "gipht.in")

% Modifications
% 2009-MAR-25: Lee, version 1
% 2009-MAR-28: Kurt, add loncenter, latcenter, tquake
% 2009-MAR-29: Kurt: add ianneal
% 2009-APR-02: Lee, revised parameter parsing to improve robustness
% 2009-MAY-10: Kurt: allow pselect == 4
% 2009-JUN-16: Kurt: allow pselect == 5
% 2009-JUL-14: Kurt: change name, add fnparin, txtoutname, objfun, fitfun, demdescfile,
% 2009-DEC-08: Kurt: allow pselect = 7
% 2010-JUN-27: Kurt: add maxcmd
% 2010-OCT-05: Kurt: add datafilename
% 2011-APR-02: Kurt: add nsaruns
% 2011-APR-18: Kurt: add parmfilename
% 2011-JUN-14: Kurt: saopt6
% 2011-JUL-03: Kurt: figopt
% 2011-JUL-23: Kurt: maxpix
% 2014-MAR-16: Kurt: surrogate

fprintf(1,'%s begin reading input control file named %s...\n',mfilename,OPT.fning);

fid = fopen(OPT.fning,'r');
if fid <= 0 
    error(sprintf ('ERROR: Cannot open input control file called %s\n',OPT.fning));
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
%     error(sprintf ('ERROR: Differing count %d keys and %d values from %s\n',nk,nv,OPT.fning));
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
      if numel(i1) > 0
          if i1 > 0 && i1 < i2
              k=k+1;
              keys{k}=strtrim(tline(1:i1-1));
              vals{k}=tline(i1+1:end);
          end
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
% fprintf(1,'%f\n',getval3(OPT,keys,vals,'xcenter','f');


% Load the values associated with keys
%  If a key isn't found it will be defaulted later
%xcenter      = getval3(OPT,keys,vals,'xcenter','f');
%if ~exist('xcenter','var'), xcenter = getval3(OPT,keys,vals,'xcenter','f'); end;
%OPT.xcenter = 9999

% xcenter      = getval3(OPT,keys,vals,'xcenter','f')   
% ycenter      = getval3(OPT,keys,vals,'ycenter','f');
% halfwidth    = getval3(OPT,keys,vals,'halfwidth','f');
% halfheight   = getval3(OPT,keys,vals,'halfheight','f');
OPT.npix         = getval3(OPT,keys,vals,'npix','f');
OPT.pselect      = getval3(OPT,keys,vals,'pselect','f');
OPT.unitv_east   = getval3(OPT,keys,vals,'unitv_east','f');
OPT.unitv_north  = getval3(OPT,keys,vals,'unitv_north','f');
OPT.unitv_up     = getval3(OPT,keys,vals,'unitv_up','f');
%tquake       = getval3(OPT,keys,vals,'tquake','f');
OPT.ithresh      = getval3(OPT,keys,vals,'ithresh','f');
OPT.ianneal      = getval3(OPT,keys,vals,'anneal','f');
OPT.pixinpatch   = getval3(OPT,keys,vals,'pixinpatch','f');
OPT.maxpix       = getval3(OPT,keys,vals,'maxpix','f');
OPT.nprocessors  = getval3(OPT,keys,vals,'nprocessors','f');
%interpcell   = getval3(OPT,keys,vals,'interpcell','f');
%mpercy       = getval3(OPT,keys,vals,'mpercy','f');
OPT.maxcmd       = getval3(OPT,keys,vals,'maxcmd','f');
OPT.nsaruns      = getval3(OPT,keys,vals,'nsaruns','f');
OPT.saopt6       = getval3(OPT,keys,vals,'saopt6','f');
OPT.figopt       = getval3(OPT,keys,vals,'figopt','b');
OPT.orbopt       = getval3(OPT,keys,vals,'orbopt','f');
OPT.surrogate    = getval3(OPT,keys,vals,'surrogate','f');
OPT.verbose      = getval3(OPT,keys,vals,'verbose','f');
% 
% % New string items 2009-JUL-14
OPT.objfun       = getval3(OPT,keys,vals,'objfun','s');
OPT.fitfun       = getval3(OPT,keys,vals,'fitfun','s');
OPT.fnparout     = getval3(OPT,keys,vals,'txtoutname','s');
OPT.demdescfile  = getval3(OPT,keys,vals,'demdescfile','s');
OPT.ilist        = getval3(OPT,keys,vals,'ilist','s');
OPT.orbfile      = getval3(OPT,keys,vals,'orbfile','s');
%cohfile      = getval3(OPT,keys,vals,'cohfile','s');
OPT.datafilename = getval3(OPT,keys,vals,'datafilename','s');
%parmfilename = getval3(OPT,keys,vals,'parmfilename','s');
OPT.printfun     = getval3(OPT,keys,vals,'printfun','s');
% pha2qlsname  = getval3(OPT,keys,vals,'pha2qlsname','s');
% phaseprefix  = getval3(OPT,keys,vals,'phaseprefix','s');
OPT.timefun     = getval3(OPT,keys,vals,'timefun','s');

% read input parameters in gin format
OPT.fnparin      = getval3(OPT,keys,vals,'fnparin','s');
% name for output parameters in gin format
OPT.fnparout     = getval3(OPT,keys,vals,'fnparout','s');
% name for output parameters in summary format
OPT.fnsumout    = getval3(OPT,keys,vals,'fnsumout','s');

% subregion to cut
OPT.cutsubregion = getval3(OPT,keys,vals,'cutsubregion','s');



% % Make sure a value is found
% %  this code can catch typos in keys
% %  or assign default values for omitted entries.
% if numel(xcenter) == 0
%     warning(sprintf ('WARNING: xcenter not yet defined in %s\n',OPT.fning));
% end
% if numel(ycenter) == 0
%     warning(sprintf ('WARNING: ycenter not yet defined in %s\n',OPT.fning));
% end
% if numel(halfwidth) == 0
% 	halfwidth = 120;	% default
% end
% if numel(halfheight) == 0
% 	halfheight = 100;	% default
% end
% if numel(npix) == 0
% 	npix = 199;	% default
% end
% if numel(pselect) == 0
% 	pselect = 1; % default to random
% end
% if numel(unitv_east) == 0
% 	unitv_east = 0.4;	% default
% end
% if numel(unitv_north) == 0
% 	unitv_north = -0.1;	% default
% end
% if numel(unitv_up) == 0
% 	unitv_up = 0.9;	% default
% end
% if numel(tquake) == 0
% 	tquake = 0;	% default
% end
% if numel(ithresh) == 0
%    ithresh = 8; % threshold circular mean deviation misfit from mean direction in DN [0 127]
% end
% if numel(maxcmd) == 0
%    maxcmd = 16; % threshold circular mean deviation from ramp in DN [0 127]
% end
% if numel(pixinpatch) == 0
%    pixinpatch = 16;   % mininum number of good pixels in a quadtree patch
% end
% if numel(maxpix) == 0
%    maxpix = 16;   % maximum number of pixels in a quadtree patch
% end
% if numel(ianneal) == 0
%    %ianneal = 1;	% default
%    ianneal = 0;	% new default 20130113
% end
% if numel(objfun) == 0
% 	objfun = str2fun('funcostrarcm'); % mininum angle, assumes zero mean, using matlab function;	% default
% end
% if numel(fitfun) == 0
%    fitfun =str2fun('funfit29');   % 
% end
% if numel(printfun) == 0
%    printfun = str2fun('printnull');   % 
% end
% if numel(fnparin) == 0
%     % try reading old key word
%     fnparin      = getval3(OPT,keys,vals,'txtfile','s');
%     if numel(fnparin) == 0
%         fnparin ='parameters.gin'; % new format
%     end
% end
% if numel(fnparout) == 0
%    fnparout=sprintf('%s_%s',runname,strrep(fnparin,'.gin','.gout'));
% end
% if numel(demdescfile) == 0
%     demdescfile = 'dem_descriptor.dat'; 
% end
% if numel(ilist) == 0
%    ilist = 'file_names.dat';
% end
if numel(OPT.nprocessors) == 0
   OPT.nprocessors = 1;
end
% if numel(interpcell) == 0
%    interpcell = 1000; % dimension of cell (in same units as xcenter and ycenter [meters or degrees]) for interpolation in gipht_step4
% end
% if numel(mpercy) == 0
%    mpercy = 0.0284; % assume ERS or ENVISAT
% end
if numel(OPT.nsaruns) == 0
    OPT.nsaruns = 1; % 1 run for SA
end
if numel(OPT.saopt6) == 0
    OPT.saopt6 = 0; % do not follow annealing with optimization
end
% if numel(figopt) == 0
%    figopt = 0; % no fancy plots
% end
% if numel(orbopt) == 0
%    orbopt = 0; % do not handle estimate orbital parameters
% end
% if numel(pha2qlsname) == 0
%     pha2qlsname ='pha2qls';   % 
% end
% if numel(phaseprefix) == 0
%     phaseprefix ='psp';   % 
% end
% if numel(surrogate) == 0
%    surrogate = 0; % use exact (rather than approximate) version of fitting function
% end
if numel(OPT.verbose) == 0
    OPT.verbose = 0; % be quiet
end
% 
% orbit file
if numel(OPT.orbfile) == 0
   %orbfile = 'master.orb';
   fprintf(1,'WARNING: No orbit file named in %s. Assuming constant look vector.\n',OPT.fning);
   % check unit vector
   OPT.unitv = [OPT.unitv_east OPT.unitv_north OPT.unitv_up];
   for i=1:3
      if abs(OPT.unitv(i)) > 1.0
         error(sprintf('ERROR: %ith component of UNITV has invalid value %f in %s\n'...
            ,i,unitv(i),OPT.fning));
      end
   end
   i=3;
%    if unitv(i) < 0
%       error(sprintf('ERROR: %ith component of UNITV has invalid value %f in %s\n'...
%          ,i,unitv(i),OPT.fning));
%    end
   unorm = abs(norm(OPT.unitv) - 1.0);
   if unorm > 1e-2
      error(sprintf('ERROR: norm of UNITV differs from unity by %f in %s\n',unorm,OPT.fning));
   end
else
   OPT.unitv = [NaN NaN NaN];
   if OPT.orbfile(end) == '/'
      fprintf(1,'Will look for individual orbit files in directory named %s\n',OPT.orbfile)
   else
      if fexist(OPT.orbfile) <= 0
         error(sprintf('Cannot find file %s as listed in %s\n',OPT.orbfile,OPT.fning));
      end
   end
end
% % coherence file
% if pselect == 3
%    if numel(cohfile) == 0
%       fprintf(1,'WARNING: No coherence file named in %s\n',OPT.fning);
%       %cohfile = 'coh.byt';
%    else
%       if fexist(cohfile) <= 0
%          error(sprintf('Cannot find file %s as listed in %s\n',cohfile,OPT.fning));
%       end
%    end
% end
% 
% 
% % Check that files exist
% Cfiles = {fnparin demdescfile ilist parmfilename};
% for i=1:numel(Cfiles);
%     if numel(Cfiles{i}) > 0
%         fprintf('%s\n',Cfiles{i});
%         if fexist(Cfiles{i}) <= 0
%             error(sprintf('Cannot find file %s as listed in %s\n',Cfiles{i},OPT.fning));
%         end
%     end
% end
% % Check that files exist
% Cfiles = {datafilename};
% for i=1:numel(Cfiles);
%     if numel(Cfiles{i}) > 0
%         fprintf('%s\n',Cfiles{i});
%         if fexist(Cfiles{i}) <= 0
%             warning(sprintf('Cannot find file %s as listed in %s\n',Cfiles{i},OPT.fning));
%         end
%     end
% end

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
cmds = {OPT.fitfun,OPT.printfun,OPT.objfun}
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

% if numel(parmfilename) == 0
%    parmfilename = 'gipht_parm.pst';
% end

% Perform sanity checks
if  ismember(OPT.pselect,[0,1,2,3,5,7,9]) == 0
    error(sprintf ('ERROR: PSELECT invalid value in %s\n',OPT.fning));
end

if  ismember(OPT.ianneal,[0:7]) == 0
    error(sprintf ('ERROR: ANNEAL invalid value in %s\n',OPT.fning));
end

% if  nsaruns <= 0  || nsaruns > 10 
%     error(sprintf ('ERROR: NSARUN invalid value in %s\n',OPT.fning));
% end

% if OPT.ithresh < 0 || OPT.ithresh > 127
%     error(sprintf('ERROR: ithresh has invalid value of %f in %s\n',OPT.ithresh,OPT.fning));
% end
% 
% if OPT.maxcmd < 0 || OPT.maxcmd > 127
%     error(sprintf('ERROR: maxcmd has invalid value of %f in %s\n',OPT.maxcmd,OPT.fning));
% end
% 
% if OPT.pixinpatch < 0 || OPT.pixinpatch > 100
%     error(sprintf('ERROR: pixinpatch has invalid value of %f in %s\n',OPT.pixinpatch,OPT.fning));
% end

% if interpcell < 0 || interpcell > 100e3
%     error(sprintf('ERROR: interpcell has invalid value of %f in %s\n',interpcell,OPT.fning));
% end

% if pselect < 3
%     if npix < 10 || npix > (2*halfwidth + 1) * (2*halfheight +1)
%         error(sprintf ('ERROR: NPIX is < 10 or > H*W in %s\n',OPT.fning));
%     end
% end

if OPT.figopt < 0 || OPT.figopt > 8
    error(sprintf ('ERROR: FIGOPT invalid value %d in %s\n',OPT.figopt,OPT.fning));
end
% % vector displacement requires interpolation
% if bitget(figopt,3) == 1 && bitget(figopt,2) == 0
%     warning(sprintf('vector displacement requires interpolation\n'));
%     figopt = bitset(figopt,2,1);
%     fprintf(1,'figopt = %4s\n',dec2bin(figopt,4));
% end

if  ~ismember(OPT.orbopt,[0,1])
    error(sprintf('ERROR: orbopt has invalid value of %f in %s\n',OPT.orbopt,OPT.fning));
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

if OPT.nprocessors > 1 
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
        current_parpool= parpool([1,OPT.nprocessors]);
    else
        error(sprintf('Request is for nprocessors = %d BUT matlab distributed tool kit not installed.\n',nprocessors));
    end
end

% if mpercy < 0.001 || mpercy > 1.000
%     error(sprintf('ERROR: mpercy has invalid value of %f in %s\n',mpercy,OPT.fning));
% end

% allowed = {'pha','psp','cln','psm'};
% if ismember(phaseprefix,allowed) ~= 1
%     allowed
%     error(sprintf('ERROR: phaseprefix %s is not on list of allowed values above.\n',phaseprefix));
% end

% subregion
if numel(OPT.cutsubregion) > 0
    % expect string as for GMT w/e/s/n for example: -R-135/-133/34/35
    cbuff = char(OPT.cutsubregion);
    cbuff = strrep(cbuff,'-R','');
    cbuff = strsplit(cbuff,'/');
    
    
else
    OPT.cutsubregion = 0; % 
end

fprintf(1,'%s completed reading input controls.\n',mfilename);
return
