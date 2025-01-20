clear
clc

%% Load files from multiple folders
main_dir = fullfile('/network/iss/cenir/analyse/irm/users/anna.skrzatek','nifti');

e_PARKGAME = exam(main_dir,'PARKGAME');
e_REMINARY = exam(main_dir,'REMINARY_\w{2}_');

e = e_PARKGAME + e_REMINARY;


%% Load files from particular folder
%main_dir = fullfile('/home/anna.skrzatek/data/','nifti');
%e = exam(main_dir,'2018_07_18_PARKGAMEII_001_NB_18_07_2018_V1');

e.addSerie('ACTIVATION$','run_ACTIVATION',1)
e.addSerie(        'RS$','run_RS'        ,1)

e.getSerie('run').addVolume('^f','f',3)

e.addSerie('t1mpr_S256_0_8iso_p2$','anat_T1',1)
e.getSerie('anat_T1').addVolume('^s.*p2.nii','s',1)


dir_func = e.getSerie('run') .toJob;
dir_anat = e.getSerie('anat').toJob(0);

par.fake = 0;
par.redo = 0;
par.verbose = 2;
par.MNI = 0; % warning : keep MNI=0 to use SPM segmentation after
par.cmd_arg = "--mask_mode='anat'";
%par.mask_mode = 'anat';
par.run = 1;


%% MEICA

par.pct = 1; % Parallel Computing Toolbox

if par.pct
    par.nrCPU = 1;
else
    par.nrCPU = 0;
end

tic
job_meica_afni(dir_func, dir_anat, par);
toc

% add new volumes
e.getSerie('run').addVolume('S\d{2}_\w+_hikts_nat\.*nii'   ,'hikts',1)
e.getSerie('run').addVolume('^S\d{2}_\w+_T1c_medn_nat\.*nii', 'T1c',1)

[ec, ei] = e.removeIncomplete;
e = ec;

e.getSerie('run').getVolume('T1c').unzip(par)
% then REVIEW !


%% Prepare CAT12 segmentation

%anat segment
fanat = e.getSerie('anat').getVolume('^s').toJob;

% Retrocompatibility for SPM:Spatial:Segment options
par.GM        = [1 1 1 1]; % warped_space_Unmodulated(wp*) / warped_space_modulated(mwp*) / native_space(p*) / native_space_dartel_import(rp*)
par.WM        = [1 1 1 1];
par.CSF       = [1 1 1 1];
par.bias      = [1 1 0] ;  % native normalize dartel     [0 1]; % bias field / bias corrected image
par.warp      = [1 1]; % warp field native->template / warp field native<-template

par.jacobian  = 0;         % write jacobian determinant in normalize space
par.doROI     = 0;
par.doSurface = 0;
par.subfolder = 0; % all results in the same subfolder

par.run     = 1;
par.display = 0;
par.pct     = 0;


%% Run CAT12 segmentation

j_segment = job_do_segmentCAT12(fanat,par);
e.getSerie('anat').addVolume('^y'  ,'y' );
e.getSerie('anat').addVolume('^ms' ,'ms');
e.getSerie('anat').addVolume('^wms','wms');


%% apply normalization on fmri MEICA volume, using SPM warp field from the segementation

% Prepare

ffunc = e.getSerie('run').getVolume('^T1c');
fy    = e.getSerie('anat').getVolume('^y');

par.display= 0;
par.run    = 1;
par.pct    = 1;


%% Run normalize

par.auto_add_obj=0;
job_apply_normalize(fy,ffunc,par);
e.getSerie('run').addVolume('^wS.*nii','wT1c',1)


e.explore

save('e','e')

