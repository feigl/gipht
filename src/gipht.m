% initialize variables
clear all;
% deal with slashes on Windows boxes
if ispc == 1
    set(0,'DefaultTextInterpreter','none');
end
close all; nf=0;format compact; 
echo off all

%splashtext = sprintf('%80s\n',help('gipht.m'));
versionnum = 2.9;
D=dir(which('gipht'));
versiondat = D.date;
versionstr = sprintf('Development version %.1f of %s',versionnum,versiondat);

splashtext = char({' '...
,' General Inversion of Phase Technique (GIPhT)'...
,' '...
,' Copyright (c) Kurt Feigl, Cliff Thurber  All rights reserved.'...
,' U.S. Patent No. 7,446,705.'...
,' '...
,' Use of this computer program [GIPhT] is governed by a license agreement'...
,' between You and the Wisconsin Alumni Research Foundation (the'...
,' "Agreement"). Unauthorized use, reproduction, or distribution of this'...
,' program or any portion of it may result in severe civil and criminal'...
,' penalties, and will be prosecuted to the maximum extent possible under'...
,' the law. By using this program, you agree to be bound by the terms of the'...
,' Agreement.'...
,' '...
,versionstr...
,' Source code:'...
,' Copyright (c) Kurt Feigl, Cliff Thurber, Lee Powell, Peter Sobol, Aaron Masters'...
,' Elena Baluyut, S. Tabrez Ali'...
,' University of Wisconsin-Madison'});
fprintf(1,'\n\n----------------   %s %s ----------\n',upper(mfilename),versionstr);
fprintf(1,'\n\n----------------   %s begins at %s ----------\n',upper(mfilename),datestr(now,31));
tstart = tic;


% Modification History below this line
% 2007-2008: Kurt
%     prototyping
% 2009-MAR-28: Lee
%     clean up for demo on Iceland subsidence: 6 pairs
% 2009-MAR-29: Kurt
%     expand to other examples
% 2009-APR-4: Lee
%     add license check and turn off the diary
% 2009-MAY-10: Kurt
%     add option to use coherence
% 2009-JUN-18:
%     Use signed char variables for all phases for speed
%     Introduce pselect = 5 for quadtree 
% 2010-JUN: Get gradients to work with pselect == 7
% 2011-JUN:
%     Fix bug in quadtree routines with midpoint of patch
%     Paramterize orbits in terms of incidence angle - correctly
%     Handle missing data
%     Handle correlated parameters
% 2011-JUL
%     Speed up step 5
%     Fix bug with gradients
% 2011-OCT
%     Clean up plots 
% 2011-OCT-11
%     for pselect == 7, use test_generalized_paretos to estimate critical
%     value of cost
% 2011-NOV-11
%     update for MATLAB R2011b
% 2012-JAN
%     gipht_step3: take max-min for parameter uncertainties
%     add bootstrap to anneal4
%     measure gradients in dimensionless strain
%     Taylor approximation
% 2012-SEP v. 2.5
%     print out derived parameters, too
%     identify bugs when step size in latitude (DL) is positive
%     Add quadtree coordinates to DST
%     Add phaseprefix to gipht.in
%     Fix vector components
% 2013-MAY v. 2.5 (for short course)
%     write_dst.m: write quad dimensions to dst_sample
%     gipht_step1.m: same as above
%     gipht_path.m: handle file separators under Windows and DOS
%     quad_tree: STILL TO DO same as above
% 2013-MAY v. 2.7 
%     improve partial derivatives
%     use comsol

% initialize paths
%giphtpath
%path('../src',path);
%path('../extern',path);

%license_check;

clockstr = clock;
runname=sprintf('%04d%02d%02d_%02d%02d%02d'...
   ,clockstr(1),clockstr(2),clockstr(3),clockstr(4),clockstr(5),round(clockstr(6)))
rundir = sprintf('%s',['x_' runname filesep 'x']);
%system(sprintf('mkdir -p x_%s',runname));
unix(sprintf('mkdir -p x_%s',runname));
runname=rundir;
diary(sprintf('%s.log',runname));

if fexist(sprintf('%s.log',runname)) ~= 1
    warning('Cannot open diary file named %s\n',sprintf('%s.log',runname));
end

% When there is no display, this returns [1 1 1 1] instead of an actual screen size. However, this relies on behavior that isn't actually specified (by the doc, for instance) to work in any particular way, so may be subject to change in the future. If you were going to use this many times, it might be wise to wrap it in a function (e.g. create an "isdisplay.m" function file), so you can easily change the implementation in the future, if needed. (This method worked as of MATLAB R2008a.) 
ss4 = get(0, 'ScreenSize');
if ss4 == [1 1 1 1 ]
%warning('off','all');
end

% Make a figure with splash message
figure;set(gca,'Visible','Off');axis([0 30 0 20]); hold on;
h=text(0.05,10,splashtext);
set(h,'FontName','Helvetica','Fontsize',10,'FontWeight','bold');
printpdf(sprintf('%s.pdf',runname));

% count the errors
nerrors = 0;

% Now run the following steps:
gipht_step1; % read in phase files, select pixels
gipht_step2; % set bounds on parameters, run simulated annealing
gipht_step3; % determine statistical uncertainties on parameter estimates
gipht_step4; % make images for quad tree
gipht_step5; % make images for entire sub-region

save_run;
telapsed = toc(tstart);
fprintf(1,'\n\n----------------   %s ended normally at %s ----------\n',upper(mfilename),datestr(now,31));
fprintf(1,    '----------------   Elapsed time: %.0f seconds\n',telapsed);

diary off
