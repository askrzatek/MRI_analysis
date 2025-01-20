%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   MULTI - ECHO    %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

main_dir = '/network/iss/cenir/analyse/irm/studies/AUDICOG/DATA/Non_chirurgicaux';
cd '/network/iss/cenir/analyse/irm/studies/AUDICOG';

e = exam(main_dir, 'AUDICOG_Suj'); % all subjects with multi-echo


%% Get files paths #matvol

% % Anat
e.addSerie('T1w$', 'anat_T1', 1 );
e.getSerie('anat').addVolume('^v_.*nii','v',1);

% Func
run_list = {'Task_nback', 'Task_localizer_audio', 'RS'};
for r = 1 : length(run_list)
    run_name = run_list{r};
    e.addSerie([run_name           '$'], ['run_' run_name], 1);
    e.addSerie([run_name '_PhysioLog$'], ['phy_' run_name], 1);
end

e.getSerie('run').addVolume('^v_.*nii$',   'v', 3);
e.getSerie('phy').addPhysio(     'dcm$', 'dcm', 1);

e.reorderSeries('name');

e.explore


%% Cluster ?

CLUSTER = 0;



%% segment cat12 #CAT12->SPM12

anat = e.gser('anat_T1').gvol('^v');

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
par.redo    = 0;
par.display = 0;

par.subfolder = 0;

par.GM        = [1 0 1 0]; % warped_space_Unmodulated (wp1*)     / warped_space_modulated (mwp1*)     / native_space (p1*)     / native_space_dartel_import (rp1*)
par.WM        = [1 0 1 0]; %                          (wp2*)     /                        (mwp2*)     /              (p2*)     /                            (rp2*)
par.CSF       = [1 0 1 0]; %                          (wp3*)     /                        (mwp3*)     /              (p3*)     /                            (rp3*)
par.TPMC      = [0 0 0 0]; %                          (wp[456]*) /                        (mwp[456]*) /              (p[456]*) /                            (rp[456]*)   This will create other probalities map (p4 p5 p6)

par.label     = [1 1 0] ;  % native (p0*)  / normalize (wp0*)  / dartel (rp0*)       This will create a label map : p0 = (1 x p1) + (3 x p2) + (1 x p3)
par.bias      = [1 1 0] ;  % native (ms*)  / normalize (wms*)  / dartel (rms*)       This will save the bias field corrected  + SANLM (global) T1
par.las       = [0 0 0] ;  % native (mis*) / normalize (wmis*) / dartel (rmis*)       This will save the bias field corrected  + SANLM (local) T1

par.warp      = [1 1];     % warp fields  : native->template (y_*) / native<-template (iy_*)

par.doSurface = 0;
par.jacobian  = 0;         % write jacobian determinant in normalize space
par.doROI     = 0;         % will compute the volume in each atlas region

job_do_segmentCAT12(anat,par);
%% AUTO ADD in jobs
% e.getSerie('anat_T1').addVolume('^p0','p0',1)
% e.getSerie('anat_T1').addVolume('^p1','p1',1)
% e.getSerie('anat_T1').addVolume('^p2','p2',1)
% e.getSerie('anat_T1').addVolume('^p3','p3',1)

%% Sort echos #MATLAB/matvol

clear par
par.run  = 1;
par.fake = 0;
par.sge  = 0;

par.redo = 0;
meinfo = job_sort_echos( e.getSerie('run') , par );


%% job_afni_proc_multi_echo #ANFI/afni_proc.py

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
par.fake = 0;
par.redo = 0;

par.seperate = 1;
par.write_nifti = 1;

par.blocks  = {'tshift', 'volreg'};

afni_prefix = char(par.blocks); 
afni_prefix = afni_prefix(:,1)';
afni_prefix = fliplr(afni_prefix);   

afni_subdir = ['afni_' afni_prefix];
par.subdir = afni_subdir;
job_afni_proc_multi_echo( meinfo, par );

e.addSerie('RS','afni','afni_RS',1)

%% do_fsl_robust_mask_epi #FSL
%% volume added in the job
%e.getSerie('run_RS').addVolume(['^' afni_prefix 'e1.nii'],[afni_prefix 'e1'],1)

fin  = e.getSerie('run').getVolume(['^' afni_prefix 'e1']);

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
par.fake  = 0;
par.redo  = 0;
par.fsl_output_format = 'NIFTI_GZ';
do_fsl_robust_mask_epi( fin, par );

% Checkpoint & unzip
par.jobname = 'unzip_and_keep__bet';
%e.getSerie('run_RS').addVolume(['^bet_Tmean_' afni_prefix
%'e1.nii$'],'bet_Tmean',1) %% volume added directly in the job

e.getSerie('run_RS').getVolume('bet_Tmean').removeEmpty().unzip_and_keep(par)


%% TEDANA #Python

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
par.fake  = 0;
par.redo  = 0;
par.pct = 0;

% cluster
par.walltime = '12:00:00';      % HH:MM:SS
par.mem      = '16G';           % ICA is very memory consuming
par.sge_nb_coeur = 2;           % I dont't know why, but 2 CPU increase the "stability" of the job on the cluster

tedana_subdir = ['tedana009a1_' afni_prefix];
job_tedana_009a1( meinfo, afni_prefix, tedana_subdir, ['bet_Tmean_' afni_prefix 'e1_mask.nii.gz'], par );


% Checkpoint & unzip
par.jobname = 'unzip_and_keep__tedana';

%e.addSerie('RS', tedana_subdir,'tedana_RS',1)
%e.getSerie('tedana_RS').addVolume('^ts_OC.nii$','ts_OC',1) %% added
%directly in the job
e.getSerie('run').getVolume('ts_OC').removeEmpty().unzip_and_keep(par)
e.getSerie('run').getVolume('Tmean').removeEmpty().unzip_and_keep(par)


%% Coregister TEDANA outputs to Anat #SPM12

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
par.redo  = 0;
par.type  = 'estimate';

src = e.getSerie('run').removeEmpty().getVolume('bet_Tmean');
oth = e.getSerie('run').removeEmpty().getVolume('(^ts_OC)|(^dn_ts_OC)');

ref = e.getSerie('run').removeEmpty().getExam.getSerie('anat_T1').getVolume('p0');

par.jobname = 'spm_coreg_epi2anat';
job_coregister(src,ref,oth,par);

e.getSerie('anat_T1').addVolume('^wp0','wp0',1) %% probably added in the job
e.getSerie('anat_T1').addVolume('^wp1','wp1',1)
e.getSerie('anat_T1').addVolume('^wp2','wp2',1)
e.getSerie('anat_T1').addVolume('^wp3','wp3',1)


%% Normalize TEDANA outputs #SPM12

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
par.redo = 0;
par.vox = [2.5 2.5 2.5]; % IMPORTANT keep original EPI voxel size
%img = e.getSerie('run').getVolume('(^ts_OC)|(^dn_ts_OC)').removeEmpty();
img = e.getSerie('run').getVolume('(^ts_OC)').removeEmpty();
img.getExam.getSerie('anat_T1').addVolume('^y','y',1)
y   = img.getExam.getSerie('anat_T1').getVolume('^y');
par.jobname = 'spm_normalize_epi';
job_apply_normalize(y,img,par);

% Nomalize Tmean, used later for PhysIO Noise ROI
img = e.getSerie('run').removeEmpty().getVolume('bet_Tmean_vte1$');
y   = img.getExam.getSerie('anat_T1').getVolume('^y');
par.jobname = 'spm_normalize_meanepi';
job_apply_normalize(y,img,par);


%%  Smooth TEDANA outputs #SPM12

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
par.redo = 0;

%e.getSerie('run').addVolume('^wts_OC','wts')
img = e.getSerie('run_RS').getVolume('wts').removeEmpty();

par.smooth   = [5 5 5];
par.prefix   = 's5';
job_smooth(img,par);
%e.getSerie('tedana_RS').addVolume('^s5wts_OC','s5wts') %% autoadd job

par.smooth   = [8 8 8];
par.prefix   = 's8';
job_smooth(img,par);
%e.getSerie('tedana_RS').addVolume('^s8wts_OC','s8wts') %% autoadd job


%% coregister WM & CSF on functionnal (using the warped mean) #SPM12
% This will be used for TAPAS:PhysIO

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
ref = e.getSerie('run');
ref = ref(:,1).getVolume('wbet_Tmean');
src = e.getSerie('anat_T1').getVolume('^wp2');
oth = e.getSerie('anat_T1').getVolume('^wp3');
par.type = 'estimate_and_write';
par.jobname = 'spm_coreg_WMCSF2wEPI';

job_coregister(src,ref,oth,par);


%% rp afni2spm #matlab/matvol

% input
%e.getSerie('run').addRP('rp','rp_afni',1)
dfile = e.getSerie('run').getRP('rp_afni').removeEmpty();

% output
output_dir = e.getSerie('run').getPath();

% go
job_rp_afni2spm(dfile, output_dir);
rp = fullfile(e.getSerie('run').getPath(),'rp_spm.txt');

%% extract physio from special dicom

% https://github.com/CMRR-C2P/MB
addpath('/network/iss/cenir/analyse/irm/studies/AUDICOG/local_toolboxs/MB')

e.getSerie('phy').getPhysio('dcm').extract()

e.getSerie('phy').getPhysio('phy').check() % takes a bit of time, use it once to verify your data


%% PhysIO nuisance regressor generation #matlab/TAPAS-PhysIO
%% Prepare files

info = e.getSerie('phy').getPhysio('info').removeEmpty;
puls = e.getSerie('phy').getPhysio('puls').removeEmpty;
resp = e.getSerie('phy').getPhysio('resp').removeEmpty;
run  = e.getSerie('run').removeEmpty;

%--------------------------------------------------------------------------
% in AFNI, outside the mask is NaN
% in SPM, outise the mask is 0
% so convert NaN to 0

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem'; 
else
    par.run = 1;
    par.sge = 0;
end
job_afni_remove_nan( run.getVolume('^wts'), par );
%--------------------------------------------------------------------------
%e.getSerie('tedana_RS').addVolume('^nwts_OC','nwts_OC',1) %% auto add in
%job
volume = run.getVolume('^nwts_OC').removeEmpty;

outdir = volume.getDir();

%rp = e.getSerie('tedana_RS').getRP('rp_spm');
%run.getExam.getSerie('anat').addVolume('^rwp2','rwp2',1) %% auto add in
%job
%run.getExam.getSerie('anat').addVolume('^rwp3','rwp3',1) %% auto add in
%job

mask = run.getExam.getSerie('anat').getVolume('^rwp[23]').squeeze;


%% Prepare job

clear par
%----------------------------------------------------------------------------------------------------------------------------------------------------
% ALWAYS MANDATORY
%----------------------------------------------------------------------------------------------------------------------------------------------------

par.physio   = 1;
par.noiseROI = 1;
par.rp       = 1;

par.TR     = 1.660;
par.nSlice = 60;

par.volume = volume;
par.outdir = outdir;

%----------------------------------------------------------------------------------------------------------------------------------------------------
% Physio
%----------------------------------------------------------------------------------------------------------------------------------------------------

par.physio_Info = info;
par.physio_PULS = puls;
par.physio_RESP = resp;

par.physio_RETROICOR        = 1;
par.physio_HRV              = 1;
par.physio_RVT              = 1;
par.physio_logfiles_vendor  = 'Siemens_Tics'; % Siemens CMRR multiband sequence, only this one is coded yet
par.physio_logfiles_align_scan = 'last';         % 'last' / 'first'
% Determines which scan shall be aligned to which part of the logfile.
% Typically, aligning the last scan to the end of the logfile is beneficial, since start of logfile and scans might be shifted due to pre-scans;
par.physio_slice_to_realign    = 'middle';       % 'first' / 'middle' / 'last' / sliceNumber (integer)
% Slice to which regressors are temporally aligned. Typically the slice where your most important activation is expected.


%----------------------------------------------------------------------------------------------------------------------------------------------------
% noiseROI
%----------------------------------------------------------------------------------------------------------------------------------------------------

par.noiseROI_mask   = mask;
par.noiseROI_volume = volume;

par.noiseROI_thresholds   = [0.95 0.70];     % keep voxels with tissu probabilty >= 95%
par.noiseROI_n_voxel_crop = [2 1];           % crop n voxels in each direction, to avoid partial volume
par.noiseROI_n_components = 10;              % keep n PCA componenets


%----------------------------------------------------------------------------------------------------------------------------------------------------
% Realignment Parameters
%----------------------------------------------------------------------------------------------------------------------------------------------------

par.rp_file = rp;

par.rp_order     = 24;   % can be 6, 12, 24
% 6 = just add rp, 12 = also adds first order derivatives, 24 = also adds first + second order derivatives
par.rp_method    = 'FD'; % 'MAXVAL' / 'FD' / 'DVARS'
par.rp_threshold = 0.5;  % Threshold above which a stick regressor is created for corresponding volume of exceeding value


%----------------------------------------------------------------------------------------------------------------------------------------------------
% Other
%----------------------------------------------------------------------------------------------------------------------------------------------------
par.print_figures = 0; % 0 , 1 , 2 , 3

% classic matvol
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem'; 
else
    par.run = 1;
    par.sge = 0;
end
par.display  = 0;
par.redo     = 0;

% cluster
par.jobname  = 'spm_physio';
par.walltime = '04:00:00';
par.mem      = '4G';


job_physio_tapas( par );


%% Save examArray

save e e
