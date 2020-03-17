function write_grids(runname, i, imA, imF, imG, imH, imE ...
    ,uns,mds,urs,ucs ...
    ,omr,omx,omy,omz ...
    ,umr,umx,umy,umz ...
    ,demx,demy ...
    ,demgrd, idatatype1)
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
% 20200304  - write description and FLIP before writing grid
INFO=grdinfo3(demgrd);
INFO.command = sprintf('Written by GIPhT run %s %s',runname, date);


switch idatatype1
    case 0  % wrapped phase
        INFO.title = sprintf('wrapped phase'); 
        INFO.zname = 'radians'
    case -1 % east component of gradient
        INFO.title = sprintf('east component of range gradient'); 
        INFO.zname = 'dimless'
     case 2   % range change in meters after unwrapping
        INFO.title = sprintf('range change'); 
        INFO.zname = 'mm'      
    otherwise
        error(sprintf('unknown idatatype1 %d\n',idatatype));
end

INFO.description = 'observed';
fn = sprintf('%s_%03d_OBSV.grd', runname, i);
grdwrite3(demx,demy,flipud(imA),fn,INFO);

INFO.description = 'modeled';
fn = sprintf('%s_%03d_MODL.grd', runname, i);
grdwrite3(demx,demy,flipud(imF),fn,INFO);

INFO.description = 'residual';
fn = sprintf('%s_%03d_RESD.grd', runname, i);
grdwrite3(demx,demy,flipud(imG),fn,INFO);

INFO.description = 'deviation';
fn = sprintf('%s_%03d_COST.grd', runname, i);
grdwrite3(demx,demy,flipud(imH),fn,INFO);

INFO.description = 'quadtree';
fn = sprintf('%s_%03d_QUAD.grd', runname, i);
grdwrite3(demx,demy,flipud(imE),fn,INFO);

% change info
INFO.title = sprintf('range change after unwrapping with model'); 
INFO.zname = 'mm';

% unwrapped observed values
INFO.description = 'observed';
fn = sprintf('%s_%03d_UOBS.grd', runname, i);
grdwrite3(demx,demy,flipud(uns),fn,INFO);

% unwrapped modeled values 
INFO.description = 'modeled';
fn = sprintf('%s_%03d_UMOD.grd', runname, i);
grdwrite3(demx,demy,flipud(mds),fn);

% unwrapped residual values
INFO.description = 'residual';
fn = sprintf('%s_%03d_URES.grd', runname, i);
grdwrite3(demx,demy,flipud(urs),fn);

% unwrapped cost values 
INFO.description = 'deviation';
fn = sprintf('%s_%03d_UDEV.grd', runname, i);
grdwrite3(demx,demy,flipud(ucs),fn);

% change info
INFO.title = sprintf('displacement after unwrapping with model'); 
INFO.description = 'observed';

% unwrapped displacement components
fn = sprintf('%s_%03d_UOMR.grd', runname, i);
grdwrite3(demx,demy,flipud(omr),fn);

% unwrapped displacement components 
fn = sprintf('%s_%03d_UOMX.grd', runname, i);
grdwrite3(demx,demy, flipud(omx), fn);

% unwrapped displacement components 
fn = sprintf('%s_%03d_UOMY.grd', runname, i);
grdwrite3(demx,demy,omy,fn);

% unwrapped displacement components 
fn = sprintf('%s_%03d_UOMZ.grd', runname, i);
grdwrite3(demx,demy,flipud(omz),fn);

INFO.description = 'modeled';
% unwrapped MODELED displacement components 
fn = sprintf('%s_%03d_UUMR.grd', runname, i);
grdwrite3(demx,demy,flipud(umr),fn);

% unwrapped MODELED displacement components
fn = sprintf('%s_%03d_UUMX.grd', runname, i);
grdwrite3(demx,demy,flipud(umx),fn);

% unwrapped MODELED displacement components 
fn = sprintf('%s_%03d_UUMY.grd', runname, i);
grdwrite3(demx,demy,flipud(umy),fn);

% unwrapped MODELED displacement components
fn = sprintf('%s_%03d_UUMZ.grd', runname, i);
grdwrite3(demx,demy,flipud(umz),fn);

return
