 function [DST,PST,TST] =simann1(objfun,bounds,options,fitfun,DST,PST,TST)
%function [p1,F,model,energy,count]=simann1(objfun,bounds,options,fitfun,DST,PST,TST)
% function [p1,F,model,energy,count]=simann1(objfun,fitfun,options,DST,PST,TST)
% 2010-OCT-11 added TST structure
% 2011-APR-04 if no output, then copy p0 onto p1
%
%
% % function [mhat,F,model,energy,count]=simann1(objfun,bounds,options...
%     ,fitfun,xyzm,tepochs,bpest,dops,DD,unitv,xd,yd,ippix1,ifast,storage ...
%     ,p,pnames,mmpercycle,fnamedat,fnamein)
%
% function [mhat,F,model,energy,count]=simann1(objfun,bounds,options...
%     ,fitfun,xyzm,tepochs,bpest,dops,DD,unitv,xd,yd,ippix1,ifast,storage ...
%     ,p,pnames,mmpercycle,fnamedat,fnamein)
% function [mhat,F,model,energy,count]=simann1(objfun,bounds,options...
%Simulated annealing algorithm that tries to find a minimum to the function 'FUN'.
%
% Use the SIMANN.f90
%
%   !   this implementation of simulated annealing was used in "global optimizatio
%   !   of statistical functions with simulated annealing," goffe, ferrier and
%   !   rogers, journal of econometrics, vol. 60, no. 1/2, jan./feb. 1994, pp.
%   !   65-100. briefly, we found it competitive, if not superior, to multiple
%   !   restarts of conventional optimization routines for difficult optimization
%   !   problems.
%
%   !   for more information on this routine, contact its author:
%   !   bill goffe, bgoffe@whale.st.usm.edu
%
%   ! this version in fortran 90 has been prepared by alan miller.
%   ! it is compatible with lahey's elf90 compiler.
%   ! n.b. the 3 last arguments have been removed from subroutine sa.   these
%   !      were work arrays and are now internal to the routine.
%   ! e-mail:  amiller @ bigpond.net.au
%   ! url   :  http://users.bigpond.net.au/amiller
%
%   ! latest revision of fortran 90 version - 3 august 1997


%INPUTS:
%
%'FUN' specifies the objective function.  This function should accept an input model vector
%as the first input a scalar cost. Additional context-specific argumentscan be passed to the
%objective function by passing them to 'ANNEAL' after 'OPTIONS'.
%
%'bounds' specifies the upper and lower limits that each model parameter can take on.  This matrix
%must have as many rows as model parameters and two columns.
%
%'OPTIONS' specifies a number of annealing options (empty matrix or zeros for defaults):
%
%     OPTIONS(1) 
%                == 1 inverse problem with simulated annealing
%                == 2 forward problem with partials 
%
%     OPTIONS(2) = ipair  (Pair number) if zero, then default to
%     small data set
%
%OUTPUTS:
%
%'mhat' is the best model found.
%
%'F' is the cost associated with the best model.
%
%'model' is a matrix containing the bestmodels after each sweep.
%
%'energy' is a vector containing the costs corresponding to the models in 'model'.
%
%'count' is the total number of models checked.
%
% INPUTS
%
%   p       == parameter vector (1 column)
%   xyzm   == easting, northing, up in m [3 rows]
%   tepochs == unique times in decimal years [1 row]
%   bpest  == orbital separation w.r.t. virutal orbit in meters [1 row]
%   dops    == doppler separation w.r.t. virtual orbit in PRF [1 row]
%   DD      == differencing matrix of 1s and 0s
%   unitv   == unit vector pointing from ground to sat: E, N, U (dimless) [3 elements]
%   ippix1  == array of pixel pointers [np x 1]
%
%
%
% get dimensions

%mparam = numel(PST);
mparam = PST.mparam

write_pst(PST,'PST.IN');
%mparam = NaN;
ndata = numel(DST);
ierr1 = 0;
ierr2 = 0;

%TST = NaN;

% ierr1 = write_fitfundat(fnamedat,fitfun,xyzm,tepochs,bpest,dops,DD,unitv...
%     ,xd,yd,ippix1,p,pnames,mmpercycle);



% ierr2 = write_fitfunin(fnamein,fitfun,mparams,p,pnames,bounds);

% if ierr1 == 0 && ierr2 == 0

%srcname = 'main.f90';
%exeext  = mexext
%exename = strrep(srcname,'.f90',sprintf('.%s',exeext(4:end)))
%exename = './a.out'
exename = './run_simman.sh'
if fexist(exename) == 1
  cmd1=exename
else
  error(sprintf('Could not find executable named %s\n',exename));
end

%if options(1) == 1
%   fnamedat = 'DST.DAT'; 
%elseif options(1) == 2
%  if options(2) > 0
%    i = options(2)
%    fnamedat = sprintf('DST_%03d.DAT',i);
%  else
%    fnamedat = 'DST.DAT';
%  end
%else 
%  error(sprintf('Unknown option %d\n'));
%end
fnamedat = 'DST.DAT'; 

% name of output file
fnameout = strrep(fnamedat,'DAT','OUT');

ierr = write_dst(DST,char(PST.fitfun),fnamedat); 

fprintf(1,'Using DST file named %s\n',fnamedat);

if fexist(fnameout) == 1
    warning('Deleting file named DST.OUT.');
    [status, result] = unix(sprintf('/bin/rm -v %s\n',fnameout))
end
if fexist('PST.OUT') == 1
    warning('Deleting file named PST.OUT.\n');
    [status, result] = unix(sprintf('/bin/rm -v %s\n','PST.OUT'))
end


% here is what the command line should look like
cmd2 = sprintf('%s PST.IN %d %d %s PST.OUT\n',fnamedat,options(1),options(2),fnameout);
fprintf(1,'Getting ready to call SIMANN\n');
cmd3=sprintf('%s %s',cmd1,cmd2);
fprintf(1,'Starting with command line:\n%s\n',cmd3);
[unixstat,unixout] = unix(cmd3);

if unixstat == 0
   fprintf(1,'successful.\n');
   unixout
   %p1 = PST.p1; % copy final estimate
end

% read final estimate
if fexist('PST.OUT') == 1
    PST = read_pst('PST.OUT');
else
    warning('Could not find file named PST.OUT. Substituting file named %s\n','PST.IN');
    [status, result] = unix(sprintf('/bin/cp -v %s %s\n','PST.IN','PST.OUT'));
    PST = read_pst('PST.OUT');
    PST.p1 = PST.p0; % copy initial estimate onto final estimate
end

if fexist(fnameout) ~= 1
    warning('Could not find file named %s . Substituting file named %s\n',fnamedat,fnameout);
    [status, result] = unix(sprintf('/bin/cp -v %s %s\n',fnamedat,fnameout));
end

DST = read_dst(fnameout);

% F=NaN;
% model=NaN;
% energy=NaN;
% count=NaN;
return

