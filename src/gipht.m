% TITLE: General Inversion of Phase Technique (GIPhT)
% 
% OVERVIEW: Synthetic aperture radar (SAR) is an active remote sensing
% technique used for measuring geophysical activity on the Earth?s surface.
% It records microwaves transmitted by a sensor (usually aboard a
% satellite) and reflected by features on the Earth?s surface (usually on
% land). The reflected signal contains information in the form of amplitude
% and phase data, and requires sophisticated post-processing. A technique
% known as interferometric SAR (InSAR) measures the difference in phase
% between two images of the same area, which can be used to measure motion
% and deformation of the ground. In most applications, the interferogram
% must be ?unwrapped? before it can be interpreted. The unwrapped
% interferogram may be used to monitor geophysical changes on the Earth?s
% surface associated with earthquakes, volcanoes, landslides or glaciers,
% or with the withdrawal of oil, gas, water or minerals by extractive
% industries. Unwrapping requires considerable computational power and
% time, and may lead to significant mistakes in the unwrapped interferogram
% and thus in its intepretation. UW-Madison researchers have developed an
% algorithm for interpreting an interferogram without the need for
% unwrapping. To do so, the invention interprets the interferogram by
% estimating parameters in a quantitative model directly from the wrapped
% phase data. Alternative unwrapping algorithms have been developed, but
% these can provide inadequate results in areas where the phase data are
% imperfect, leading to errors in the unwrapped phase values. Likewise,
% these algorithms rarely, if ever, provide uncertainty estimates, limiting
% attempts to weight the data in statistical analysis. Implementation of
% the invention would reduce the time and resources necessary for advanced
% interpretation of InSAR data products, and would provide a more accurate
% result that includes an assessment of the uncertainties of the parameter
% estimates.
% 
% APPLICATIONS: InSAR for monitoring hazardous natural phenomena, e.g.,
% landslides InSAR for monitoring subsidence due to extraction, e.g., oil,
% gas, water KEY
% 
% BENEFITS: Validated on real (noisy) data in a peer-reviewed publication
% Provides a more direct path to a quantitative interpretation than
% existing methods Provides a realistic assessment of uncertainty, unlike
% existing methods
% 
% PUBLICATIONS: Feigl, K. L., and C. H. Thurber (2009) A method for
% modelling radar interferograms without phase unwrapping: application to
% the M 5 Fawnskin, California earthquake of 1992 December 4 Geophys. J.
% Int., 176, 491-504. http://dx.doi.org/10.1111/j.1365-246X.2008.03881.x
% Interferometric analysis of synthetic aperture radar images (InSAR)
% measures the phase shifts between two images acquired at two distinct
% times. These ambiguous 'wrapped' phase values range from -1/2 to +1/2
% cycles. The standard approach interprets the phase values in terms of the
% change in distance between the ground and the radar instrument by
% resolving the integer ambiguities in a process known as 'unwrapping'. To
% avoid unwrapping, we have developed, validated and applied a new method
% for modelling the wrapped phase data directly. The method defines a cost
% function in terms of wrapped phase to measure the misfit between the
% observed and modelled values of phase. By minimizing the cost function
% with a simulated annealing algorithm, the method estimates parameters in
% a non-linear model. Since the wrapped phase residuals are compatible with
% a von Mises distribution, several parametric statistical tests can be
% used to evaluate the fit of the model to the data. The method, named
% General Inversion for Phase Technique (GIPhT), can handle noisy, wrapped
% phase data. Applying GIPhT to two interferograms in the area of Fawnskin,
% California, we estimate a set of model parameters describing a magnitude
% 5 aftershock of the 1992 Landers earthquake. The resulting simulation
% fits the data well. The phase final residuals have a circular mean
% deviation less than 0.15 cycles per datum. Sampling the final residuals,
% we find the circular standard deviation of a phase measurement to be
% approximately 0.2 cycle, corresponding to 6 mm in range.
% 
% LICENSING:
% 
% Copyright (c) 2009 Kurt Feigl, Cliff Thurber All rights reserved. 
%
% U.S. Patent No. 7,446,705.
% 
% Copyright (c) Kurt Feigl, Cliff Thurber, Lee Powell, Peter Sobol, Aaron
% Masters, S. Tabrez Ali, Elena C. Baluyut, University of Wisconsin-Madison
% 
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU LESSER GENERAL PUBLIC LICENSE as published by
% the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU LESSER GENERAL PUBLIC LICENSE
% along with this program. If not, see
% http://www.gnu.org/licenses/lgpl.html.

% initialize variables
clear all;
% deal with slashes on Windows boxes
if ispc == 1
    set(0,'DefaultTextInterpreter','none');
end
close all; nf=0;format compact; 
echo off all

%splashtext = sprintf('%80s\n',help('gipht.m'));

fprintf(1,'\n\nGeneral Inversion of Phase Technique (GIPhT)\n\n');
versionnum = 2.92;
D=dir(which('gipht'));
versiondat = D.date;
versionstr = sprintf('GIPhT Development version %.1f of %s'...
    ,versionnum,versiondat);

help gipht_splash
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
% 2015-JUN v. 2.9.1
%     first version in public repository on GitHub with 
%      Lesser GNU Public License
% 2015-JUN v. 2.9.2
%     handle Quadtree data in 14-column format from JPL

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

% % Make a figure with splash message
% figure;set(gca,'Visible','Off');axis([0 30 0 20]); hold on;
% h=text(0.05,10,char(splashtext));
% set(h,'FontName','Helvetica','Fontsize',10,'FontWeight','bold');
% printpdf(sprintf('%s.pdf',runname));

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
