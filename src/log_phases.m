function status = log_phases(runname, i, imA, imF, imG, imH, imE ...
    ,uns,mds,urs,ucs...
    ,omr,omx,omy,omz...
    ,umr,umx,umy,umz)
% status = log_phases(runname, i, imA, imF, imG, imH, imE)
% 24mar09lap added writing phases for this runname
% Write_pha  imA = OBSV  imF = MODL  imG = RESD  imH = COST
% Input data are scaled up by 256 for saving as int8 [-128,127]
% 2010-NOV-06 Kurt: add imE for quad-tree resampled phase
% 2011-JUN-17 - Confirm that phases are in cycles
% 2012-MAY-10 - write unwrapped quantities in mm
% 2012-JUN-26 - write unwrapped vector quantities in mm
% 2012-JUN-27 - write unwrapped modeled vector quantities in mm

scale=256;

if strcmp(getenv('NOPHASES'),'TRUE') == 1
    disp('Env variable NOPHASES is set to TRUE, skipping log_phases');
    return
end

fn = sprintf('%s_%03d_OBSV.pha', runname, i);

status = write_pha(fn, scale*imA);
if status ~= 0
	fprintf(1,'Trouble writing %s',fn);
	return
end

fn = sprintf('%s_%03d_MODL.pha', runname, i);
status = write_pha(fn, scale*imF);
if status ~= 0
	fprintf(1,'Trouble writing %s',fn);
	return
end

fn = sprintf('%s_%03d_RESD.pha', runname, i);
status = write_pha(fn, scale*imG);
if status ~= 0
	fprintf(1,'Trouble writing %s',fn);
	return
end

% Write costs as signed bytes even though they are always positive.
fn = sprintf('%s_%03d_COST.pha', runname, i);
status = write_pha(fn, scale*imH);
if status ~= 0
	fprintf(1,'Trouble writing %s',fn);
	return
end

fn = sprintf('%s_%03d_QUAD.pha', runname, i);
status = write_pha(fn, scale*imE);
if status ~= 0
	fprintf(1,'Trouble writing %s',fn);
	return
end

% unwrapped observed values in MILLIMETERS
fn = sprintf('%s_%03d_UOBS.i2', runname, i);
status = write_i2(fn, uns);
if status ~= 0
	fprintf(1,'Trouble writing %s',fn);
	return
end

% unwrapped modeled values in MILLIMETERS
fn = sprintf('%s_%03d_UMOD.i2', runname, i);
status = write_i2(fn, mds);
if status ~= 0
	fprintf(1,'Trouble writing %s',fn);
	return
end

% unwrapped residual values in MILLIMETERS
fn = sprintf('%s_%03d_URES.i2', runname, i);
status = write_i2(fn, urs);
if status ~= 0
	fprintf(1,'Trouble writing %s',fn);
	return
end

% unwrapped cost values in MILLIMETERS
fn = sprintf('%s_%03d_UDEV.i2', runname, i);
status = write_i2(fn, ucs);
if status ~= 0
	fprintf(1,'Trouble writing %s',fn);
	return
end

% 2012-JUN-26
% unwrapped displacement components in MILLIMETERS
fn = sprintf('%s_%03d_UOMR.i2', runname, i);
status = write_i2(fn, omr);
if status ~= 0
	fprintf(1,'Trouble writing %s',fn);
	return
end
% unwrapped displacement components in MILLIMETERS
fn = sprintf('%s_%03d_UOMX.i2', runname, i);
status = write_i2(fn, omx);
if status ~= 0
	fprintf(1,'Trouble writing %s',fn);
	return
end
% unwrapped displacement components in MILLIMETERS
fn = sprintf('%s_%03d_UOMY.i2', runname, i);
status = write_i2(fn, omy);
if status ~= 0
	fprintf(1,'Trouble writing %s',fn);
	return
end
% unwrapped displacement components in MILLIMETERS
fn = sprintf('%s_%03d_UOMZ.i2', runname, i);
status = write_i2(fn, omz);
if status ~= 0
	fprintf(1,'Trouble writing %s',fn);
	return
end

% 2012-JUN-27
% unwrapped MODELED displacement components in MILLIMETERS
fn = sprintf('%s_%03d_UUMR.i2', runname, i);
status = write_i2(fn, umr);
if status ~= 0
	fprintf(1,'Trouble writing %s',fn);
	return
end
% unwrapped MODELED displacement components in MILLIMETERS
fn = sprintf('%s_%03d_UUMX.i2', runname, i);
status = write_i2(fn, umx);
if status ~= 0
	fprintf(1,'Trouble writing %s',fn);
	return
end
% unwrapped MODELED displacement components in MILLIMETERS
fn = sprintf('%s_%03d_UUMY.i2', runname, i);
status = write_i2(fn, umy);
if status ~= 0
	fprintf(1,'Trouble writing %s',fn);
	return
end
% unwrapped MODELED displacement components in MILLIMETERS
fn = sprintf('%s_%03d_UUMZ.i2', runname, i);
status = write_i2(fn, umz);
if status ~= 0
	fprintf(1,'Trouble writing %s',fn);
	return
end

return
