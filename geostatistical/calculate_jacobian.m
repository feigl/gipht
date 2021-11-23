function [nJ,mJ] = calculate_jacobian(igflib,folderName_base,pointing,GFUN_fname,VOXX_fname,NODE_fname,META_fname)
%function ierr = calculate_jacobian(igflib,folderName_base,pointing,GFUN_fname,VOXX_fname,NODE_fname,META_fname)
% calculate partial derivative of LOS displacement with respect to unit volumetric strain
% 2021/09/11 Kurt Feigl and Sui "Jay" Tung
% initialize
nJ=0;
mJ=0;
        
switch igflib
    case {3,4,5}
        l0 = 900; % element size in meters
        w0 = 900;  % original dimension in Y [m]
        h0 = 900;  % original dimension in Z [m]
        % each square expands by this change in dimensions in meters
        dl = 0.5;
        % original volume
        V0 = l0 * w0 * h0; % cubic meters
        % edge of square pyramid base
        edge0=sqrt(l0^2*2)
        % volume of right rectangular pyramid V = lwh/3
        % total volume change for 6 pyramids
        DVG = 6 * (dl * edge0 * edge0)/3.
        % second volume
        V1 = V0 + DVG

        % get list of all files
        folderName_x=strcat(folderName_base,'x');
        folderName_y=strcat(folderName_base,'y');
        folderName_z=strcat(folderName_base,'z');
        
        % sort list of files in "natural" order, even without leading zeros
        DIR_z=natsortfiles(dir(strcat(folderName_z,filesep,filesep,'gfdata_nodepair',filesep,'coord_displacement_*')));nfiles_z = numel(DIR_z);
        DIR_x=natsortfiles(dir(strcat(folderName_x,filesep,filesep,'gfdata_nodepair',filesep,'coord_displacement_*')));nfiles_x = numel(DIR_x);
        DIR_y=natsortfiles(dir(strcat(folderName_y,filesep,filesep,'gfdata_nodepair',filesep,'coord_displacement_*')));nfiles_y = numel(DIR_y);
        
        % check number of files
        if nfiles_x == nfiles_z && nfiles_y == nfiles_z
            nvoxels = nfiles_z
        else
            nfiles_x
            nfiles_y
            nfiles_z
            warning('number of files do not agree. Using minimum/n');
            nvoxels = min([nfiles_x,nfiles_y,nfiles_z])
        end
        
        % Define coordinates of VOXELS in reservoir
        % Z is with respect to ground surface positive UPWARDS
        % 2021/06/22 note extra "z" in name of folder
        %VOXX = readtable('/System/Volumes/Data/mnt/t31/jay/calculixjob/COSO_05282022_flat_HOM_500_v3/z/planar_set_id_z_xyz_pool.3dgf','filetype','text');
        VOXX = readtable(strcat(folderName_base,filesep,'z',filesep,'planar_set_id_z_xyz_pool.3dgf'),'filetype','text');
        nelements = numel(VOXX)
        VOXX.Properties.VariableNames = {'nodeid','X','Y','Z'}; % in meters
        % Truncate to those voxels within reservoir domain
        VOXX=VOXX(1:nvoxels,:);
        
        % These coordinates are with respect to origin at this location
        MESH = load(strcat(folderName_base,filesep,'mpf_fem_center_fem_utm.slave1'),'filetype','text')
    otherwise
        error('unknown igflib')
end

%% Calculate displacement field by applying weighting to Green's functions
% In the folder, the file corresponds to the displacement signature of unit dilatation of each
% lateral face for the ith voxel. Within the file, the first 3 columns are
% the top surface node location (w.r.t the fem model center documented in
% mpf_fem_center_fem_utm.slave1) and the last 3 columns are dx, dy, and dz
% respectively.

% Total LOS displacment
ULOStotal = 0.;
icol = 0;

for i=1:nvoxels
    disp([i,nvoxels])
    % get displacement at every pixel for a model with only one source voxel
    fileName_z = DIR_z(i).name;
    if contains(fileName_z,sprintf('%d',i))
        TGF1_z = readtable(strcat(folderName_z,filesep,'gfdata_nodepair',filesep,fileName_z),'FileType','text');
        TGF1_z.Properties.VariableNames={'X','Y','Z','Ux','Uy','Uz'};
    else
        i
        error('mismatched file name in Z');
    end
    
    fileName_x = DIR_x(i).name;
    if contains(fileName_x,sprintf('%d',i))
        TGF1_x = readtable(strcat(folderName_x,filesep,'gfdata_nodepair',filesep,fileName_x),'FileType','text');
        TGF1_x.Properties.VariableNames={'X','Y','Z','Ux','Uy','Uz'};
    else
        i
        error('mismatched file name in X');
    end
    
    fileName_y = DIR_y(i).name;
    if contains(fileName_y,sprintf('%d',i))
        TGF1_y = readtable(strcat(folderName_y,filesep,'gfdata_nodepair',filesep,fileName_y),'FileType','text');
        TGF1_y.Properties.VariableNames={'X','Y','Z','Ux','Uy','Uz'};
    else
        i
        error('mismatched file name in Y');
    end
    
    % column vector with one row for each observation point (pixel in interferogram)
    icol = icol+1;
    
    % start building columns
    if icol == 1
        Node1X = TGF1_z.X;
        Node1Y = TGF1_z.Y;
        Node1Z = TGF1_z.Z;
    else
        % check that locations of nodes match
        if   numel(find(abs(TGF1_z.X - Node1X) > eps)) > 0 ...
                || numel(find(abs(TGF1_z.Y - Node1Y) > eps)) > 0 ...
                || numel(find(abs(TGF1_z.Z - Node1Z) > eps)) > 0
            i
            TGF1_z.X(i)
            Node1X
            TGF1_z.Y(i)
            Node1Y
            TGF1_z.Z(i)
            Node1Z
            error('Coordinates of nodes do not match.');
        end
    end
    
    % displacement is sum of each face
    displacement_X=(TGF1_x.Ux+TGF1_y.Ux+TGF1_z.Ux);
    displacement_Y=(TGF1_x.Uy+TGF1_y.Uy+TGF1_z.Uy);
    displacement_Z=(TGF1_x.Uz+TGF1_y.Uz+TGF1_z.Uz);
    
    % LOS displacement for volume change of 1 cubic meter
    % 2021/09/04 Comment above is INCORRECT
    % 2021/09/04 Following code is LOS displacement [m]
    % 2021/09/04 for one unit of volumetric strain
    % Dimensions of GF are meters/strain
    % TODO consider using microstrain for conditioning?
    % LOS is reckoned positive toward the satellite
    % LOS is negative for subsidence
    ULOS =+1*((displacement_X/(DVG/V0)) * pointing(1) ...
        +     (displacement_Y/(DVG/V0)) * pointing(2) ...
        +     (displacement_Z/(DVG/V0)) * pointing(3));
    
    % start building table for Jacobian matrix of partial derivative of
    % observable quantity with respect to parameters
    % each row corresponds to the location of one pixel
    % each colum corresponds to a parameter
    if icol == 1
        GFUN = table(colvec(ULOS));
        GFUN.Properties.VariableNames = {fileName_z};
    else
        GFUN = [GFUN, table(colvec(ULOS))];
        GFUN.Properties.VariableNames{icol} = fileName_z;
    end
    
    % sum displacement field over all source voxels
    % assuming linear superposition
    ULOStotal = ULOStotal + ULOS;
end

% check results
[nJ,mJ] = size(ULOStotal);

%% build a table of metadata
META.pointing = pointing;
META.basefoldername = folderName_base;
META.nvoxels=nvoxels;
META.l0=l0;
META.w0=w0;
META.h0=h0;
META.meshx=MESH(1);
META.meshy=MESH(2);
META.meshz=MESH(3);
META.DV = DVG;
META.V0 = V0;
META.V1 = V1;

% write four separate tables
%     writetable(TGF,TGF_fname);
%     writetable(VOXX,TVOX_fname);
%     writetable(TGF1_z,NODE_fname);
%     writetable(TMETA,TMETA_fname);
% write four separate tables
save(GFUN_fname,'GFUN','-v7.3');
save(VOXX_fname,'VOXX','-v7.3');
%save(NODE_fname,'TGF1_z','-v7.3');
NODE = TGF1_z;
clear TGF1_z;  % change name
save(NODE_fname,'NODE','-v7.3');
save(META_fname,'META','-v7.3');
return
end


