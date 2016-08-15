function write_grids(runname, i, imA, imF, imG, imH, imE ...
    ,uns,mds,urs,ucs...
    ,omr,omx,omy,omz...
    ,umr,umx,umy,umz...
    ,demx,demy)
% function write_grids(runname, i, imA, imF, imG, imH, imE ...
%     ,uns,mds,urs,ucs...
%     ,omr,omx,omy,omz...
%     ,umr,umx,umy,umz...
%     ,demx,demy)
% status = log_phases(runname, i, imA, imF, imG, imH, imE)
% 24mar09lap added writing phases for this runname
% Write_pha  imA = OBSV  imF = MODL  imG = RESD  imH = COST
% Input data are scaled up by 256 for saving as int8 [-128,127]
% 2010-NOV-06 Kurt: add imE for quad-tree resampled phase
% 2011-JUN-17 - Confirm that phases are in cycles
% 2012-MAY-10 - write unwrapped quantities in mm
% 2012-JUN-26 - write unwrapped vector quantities in mm
% 2012-JUN-27 - write unwrapped modeled vector quantities in mm


fn = sprintf('%s_%03d_OBSV.grd', runname, i);
grdwrite3(demx,demy,imA,fn);

fn = sprintf('%s_%03d_MODL.grd', runname, i);
grdwrite3(demx,demy,imF,fn);

fn = sprintf('%s_%03d_RESD.grd', runname, i);
grdwrite3(demx,demy,imG,fn);

% Write costs as signed bytes even though they are always positive.
fn = sprintf('%s_%03d_COST.grd', runname, i);
grdwrite3(demx,demy,imH,fn);

fn = sprintf('%s_%03d_QUAD.grd', runname, i);
grdwrite3(demx,demy,imE,fn);

% unwrapped observed values 
fn = sprintf('%s_%03d_UOBS.grd', runname, i);
grdwrite3(demx,demy,uns,fn);

% unwrapped modeled values 
fn = sprintf('%s_%03d_UMOD.grd', runname, i);
grdwrite3(demx,demy,mds,fn);

% unwrapped residual values 
fn = sprintf('%s_%03d_URES.grd', runname, i);
grdwrite3(demx,demy,urs,fn);

% unwrapped cost values in
fn = sprintf('%s_%03d_UDEV.grd', runname, i);
grdwrite3(demx,demy,ucs,fn);

% unwrapped displacement components
fn = sprintf('%s_%03d_UOMR.grd', runname, i);
grdwrite3(demx,demy,omr,fn);

% unwrapped displacement components 
fn = sprintf('%s_%03d_UOMX.grd', runname, i);
grdwrite3(demx,demy, omx, fn);

% unwrapped displacement components 
fn = sprintf('%s_%03d_UOMY.grd', runname, i);
grdwrite3(demx,demy,omy,fn);

% unwrapped displacement components 
fn = sprintf('%s_%03d_UOMZ.grd', runname, i);
grdwrite3(demx,demy,omz,fn);

% unwrapped MODELED displacement components 
fn = sprintf('%s_%03d_UUMR.grd', runname, i);
grdwrite3(demx,demy, umr,fn);

% unwrapped MODELED displacement components
fn = sprintf('%s_%03d_UUMX.grd', runname, i);
grdwrite3(demx,demy,umx,fn);

% unwrapped MODELED displacement components 
fn = sprintf('%s_%03d_UUMY.grd', runname, i);
grdwrite3(demx,demy,umy,fn);

% unwrapped MODELED displacement components
fn = sprintf('%s_%03d_UUMZ.grd', runname, i);
grdwrite3(demx,demy,umz,fn);

return
