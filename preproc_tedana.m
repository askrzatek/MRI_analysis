% preproc_tedana

clc
clear all

%% Initialisation
CLUSTER = 1;
par.pct = 0;

main_dir = fullfile('/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek','nifti_test');
%main_dir = fullfile('/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test', 'firstlevel_sts_tapas_PARK_doublerun'); %,'PRISMA_REMINARY'); % pour comparer les résultats avec les VS

e_PARKGAME = exam(main_dir,'PARKGAME.*._[a,c]$');
e_REMINARY = exam(main_dir,'REMINARY_\w{2}_');

%e_PARKGAME = exam(main_dir,'PARKGAME.*._exclu$');
e = e_PARKGAME; % (3:length(e_PARKGAME)); % choose specific
%e = e_REMINARY;
e.addSerie('ACTIVATION$','run_ACTIVATION',1)
e.addSerie(        'RS$','run_RS'        ,1)

%e.getSerie('run').addVolume('^f\d{3}','f',3)
e.getSerie('run_ACTIVATION').addVolume('^v.*ACTIVATION.*nii','f',3)
e.getSerie('run_RS').addVolume('^v.*RS.*nii','f',3)
%e.addSerie('t1mpr_.*p2$','anat_T1',1)
e.addSerie('t1mpr_S256_0_8iso_p2$','anat_T1',1)
e.getSerie('anat_T1').addVolume('^v.*p2.nii','s',1)

e.reorderSeries('path');
e.unzipVolume(par);

dir_func = e.getSerie('run') .toJob;
dir_anat = e.getSerie('anat').toJob(0);


%% vol quantity control %% not useful anymore - file names changed with xnat - need to update if want to use this functionality
% % functional ACTI
% p = e.getSerie('run_ACTIVATION').getVolume('^f').print;
% p = cellstr(p);
% res =regexp(p,'f630');
% res = cellfun(@isempty,res);
% bad_volumes = char(p(res));
% 
% if ~isempty(bad_volumes)
%     disp(bad_volumes)
%     %error('wrong number of volumes')
% else
%     fprintf('all echos seem to have 630 volumes \n')
% end
% 
% % functional RS
% p = e.getSerie('run_RS').getVolume('^f').print;
% p = cellstr(p);
% res =regexp(p,'f300');
% res = cellfun(@isempty,res);
% bad_volumes = char(p(res));
% 
% if ~isempty(bad_volumes)
%     disp(bad_volumes)
%     %error('wrong number of volumes')
% else
%     fprintf('all echos seem to have 300 volumes \n')
% end

%% sort echos

clear par

if CLUSTER
    par.run      = 0;
    par.sge      = 1;
else
    par.run      = 1;
    par.sge      = 0;
end
par.fname        = 'meinfo'; % name of the .mat file that will be saved
par.redo         = 1;

par.fake         = 0;

run_list = e.getSerie('run');

[meinfo, jobs] = job_sort_echos( run_list , par )

e.gser('run').addVolume('^e','e',3);

meinfo.volume = e.getSerie('run').getVolume('^e');
meinfo.anat   = e.getSerie('anat').getVolume('^s');

%% AFNI - vtd

clear par

par.blocks   = {'despike','tshift','volreg'}; % now codded : despike, tshift, align, volreg
par.seperate    = 1;                             % each volume is treated seperatly : useful when runs have different orientations
par.execute     = 1;                             % execute afni_proc.py generated tcsh script file immidatly after the generation
par.write_nifti = 1;
if CLUSTER
    par.run     = 0;
    par.sge     = 1;
else
    par.run     = 1;
    par.sge     = 0;
end

par.jobname  = 'job_afni_proc_multi_echo';
par.subdir   = 'afni';

par.pct      = 0;

par.redo     = 0;
par.fake     = 0;
par.verbose  = 1;

tic
job_afni_proc_multi_echo(meinfo, par)
toc

%% add new vols to e

e.getSerie('run').addVolume('^vtde1.nii','vtde1',1)
e.getSerie('run').addVolume('^vtde2.nii','vtde2',1)
e.getSerie('run').addVolume('^vtde3.nii','vtde3',1)

[ec_vtd, ei_vtd ] = e.removeIncomplete; % need to check if is working as supposed
e = ec_vtd;

%% do_fsl_robust_mask_epi

fin = e.getSerie('run').gvol('^vtde1');

clear par
if CLUSTER
    par.run    = 0;
    par.sge    = 1;
else
    par.run    = 1;
    par.sge    = 0;
end

par.fake       = 0;
par.redo       = 0;

par.fsl_output_format = 'NIFTI_GZ';
[job, fmask] = do_fsl_robust_mask_epi(fin, par);
%[job, fmask] = do_fsl_robust_mask_epi({fin.gpath{:}}, par);

% copy original gz & unzip
e.gser('run').addVolume('^bet.*vtde1.nii.gz','bet',1);
e.gser('run').gvol('bet').removeEmpty.unzip_and_keep(par);
%e.gser('run').addVolume('^bet.*vtde1.nii$','bet',1);


% copy original gz mask & unzip
e.gser('run').addVolume('^bet.*vtde1_mask.nii','mask',1);
e.gser('run').gvol('mask').removeEmpty.unzip_and_keep(par);
%e.gser('run').addVolume('^bet.*vtde1_mask.nii$','mask',1);

% %% temporal mean #obsolete !!!
% 
% clear par
% par.run   = 1;
% par.fake  = 0;
% par.sge   = 0;
% par.fsl_output_format = 'NIFTI';
% 
% par.jobname = 'fslmean';
% vtde_in  = e.getSerie('run').getVolume('^vtde1').removeEmpty.toJob;
% vtde_out = addprefixtofilenames(vtde_in,'mean_');
% do_fsl_mean(vtde_in,vtde_out,par);
% 
% e.getSerie('run').addVolume('^mean_vtde1.nii','mean_vtde1',1)
% 
% 
% %% bet (brain extraction)
% 
% clear par
% par.run   = 1;
% par.fake  = 0;
% par.sge   = 0;
% par.fsl_output_format = 'NIFTI';
% 
% run_dir = e.getSerie('run').removeEmpty.toJob(0);
% 
% job = cell(0);
% for iJob = 1 : length(run_dir)
%     job{iJob,1} = sprintf('export FSLOUTPUTTYPE=NIFTI; cd %s; /network/lustre/iss02/cenir/software/irm/fsl5/bin/bet %s %s -R -m -f 0.3',...
%         run_dir{iJob},...
%         'mean_vtde1.nii',...
%         'bet_mean_vtde1.nii');
% end
% par.jobname = 'fslbet';
% do_cmd_sge(job,par);
% 
% e.getSerie('run').addVolume('^bet_mean_vtde1_mask.nii','mask_vtde1',1)
% e.getSerie('run').addVolume('^bet_mean_vtde1.nii'     , 'bet_vtde1',1)
% 

%% Tedana

% tedana.py main arguments
par.tedpca     = 'mle'; % mdl, kundu, kundu-stabilize (tedana default = mdl)
par.maxrestart = []; % (tedana default = 10)
par.maxit      = []; % (tedana default = 500)
%par.png        = 1;  % tedana will make some PNG files of the components Beta map, for quality checks

% tedana.py other arguments
par.cmd_arg = ''; % Allows you to use all addition arguments not scripted in this job_tedana.m file

% matvol classic options

par.pct      = 0; % Parallel Computing Toolbox, will execute in parallel all the subjects
par.sge      = 1;
par.redo     = 0; % overwrite previous files
par.fake     = 0; % do everything exept running
par.verbose  = 2; % 0 : print nothing, 1 : print 2 first and 2 last messages, 2 : print all
par.run      = 0;

% Cluster
if par.sge               % for ICM cluster, run the jobs in paralle
    par.jobname      = 'job_tedana';
    par.walltime     = '12:00:00';      % HH:MM:SS
    par.mem          = '16G';           % MB
    %par.sge_queu = 'normal,bigmem'; % use both
    par.sge_nb_coeur = 2;
end

%!source python_path3.6
%!source activate tedana_0.0.9a1

job_tedana_009a1( meinfo, 'vtd', 'tedana009a1_vtd', 'bet_Tmean_vtde1_mask.nii.gz', par );


%% tedana report % meinfo is in main_dir and tedana outputs in the subdir of run dir (until change - cannot use tedana_report)
par.subdir = 'tedana009a1_vtd';
%tedana_report (main_dir, par); % changement et necessite de meinfo.path

%% Add new preprocessed data to exam object : dn_ts_OC.nii for tedana denoising and ts_OC.nii for TAPAS denoising

e.addSerie('ACTIVATION','tedana009.*_vtd','tedana_ACTIVATION',1);
e.addSerie('RS','tedana009.*_vtd','tedana_RS',1);
e.getSerie('tedana').addVolume('^ts_OC','ts',1);     % tapas
%e.getSerie('tedana').addVolume('^dn_ts','dn_ts',1); % tedana

% copy original gz & unzip 
e.getSerie('tedana').getVolume('.*ts').removeEmpty.unzip_and_keep(par)

% if already unziped then
e.getSerie('tedana').addVolume('^ts_OC.nii','ts',1); % tapas

% unzip the outputs of tedana
%e.unzipVolume();
% lines below surely not needed anymore
%[ec_dn, ei_dn] = e.removeIncomplete;
%e = ei_dn;

% if all in run directory instead of tedana_vtd_mle
%e(1).getSerie('run').addVolume('^dn.*.nii$','dn',1);

%% Run CAT12 segmentation
clear par
%anat segment
fanat = e.getSerie('anat').getVolume('^s').getPath;

% Retrocompatibility for SPM:Spatial:Segment options
par.GM        = [1 0 1 1]; % warped_space_Unmodulated(wp*) / warped_space_modulated(mwp*) / native_space(p*) / native_space_dartel_import(rp*)
par.WM        = [1 0 1 1];
par.CSF       = [1 0 1 1];
par.TPMC      = [0 0 0 0];

par.bias      = [1 1 0] ;  % native normalize dartel     [0 1]; % bias field / bias corrected image
par.warp      = [1 1]; % warp field native->template / warp field native<-template
par.label     = [1 1 0];
par.las       = [0 0 0];

par.jacobian  = 0;         % write jacobian determinant in normalize space
par.doROI     = 0;
par.doSurface = 0;
par.subfolder = 0; % all results in the same subfolder

par.display = 0;
par.pct     = 0;
par.redo    = 0;

if CLUSTER
    par.run     = 0;
    par.sge     = 1;
else
    par.run     = 1;
    par.sge     = 0;
end

j_segment = job_do_segmentCAT12(fanat,par);


e.getSerie('anat').addVolume('^y'  ,'y' );
e.getSerie('anat').addVolume('^p0' ,'p0' );

%% SPM
%% Coregister TEDANA outputs to Anat

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
else
    par.run = 1;
    par.sge = 0;
end

% ACTIVATION
clear par
par.type   = 'estimate';
par.interp = 1;
par.prefix = 'r';
par.redo   = 0;
par.display= 0;
par.cost_fun = 'ncc';

tmp_exam_a = [e.gser('run_ACTIVATION').removeEmpty.exam];
ref = tmp_exam_a.gser('anat').getVolume('^p0');
e.getSerie('run_ACTIVATION').getVolume('bet').removeEmpty.unzip_and_keep;
src = e.getSerie('run_ACTIVATION').getVolume('bet');
oth = e.getSerie('tedana_ACTIVATION').removeEmpty.getVolume('.*ts'); % on peut qjouter le masque pour qu'il suive l'anat

%keep in mind this option if emergency to get the precedent header
%do_fsl_copy_header(influencer, origin,      par);
%do_fsl_copy_header(follower,   influencer,  par);


job_coregister(src,ref,oth,par);

%apply same transformation to the mask as the one for bet_Tmean
follower = e.getSerie('run_ACTIVATION').removeEmpty.getVolume('mask') .toJob;
influencer = oth. toJob;
do_fsl_copy_header(follower,   influencer,  par);
   
% RS
tmp_exam_rs = [e.getSerie('run_RS').removeEmpty.exam];
ref = tmp_exam_rs.getSerie('anat').getVolume('^p0');
e.getSerie('run_RS').removeEmpty.getVolume('bet').removeEmpty.unzip_and_keep;
src = e.getSerie('run_RS').removeEmpty.getVolume('bet');
oth = e.getSerie('tedana_RS').removeEmpty.getVolume('^ts'); %changed from .* to ^

job_coregister(src,ref,oth,par);
    
%apply same transformation to the mask as the one for bet_Tmean
follower = e.getSerie('run_RS').removeEmpty.getVolume('mask') .toJob;
influencer = oth. toJob;
do_fsl_copy_header(follower,   influencer,  par);


%% Normalize TEDANA outputs

clear par
par.redo         = 0;
if CLUSTER
    par.sge      = 1;
    par.run      = 0;
else
    par.sge      = 0;
    par.run      = 1;
end
par.auto_add_obj = 0;
par.display      = 0;
par.vox          = [2.5 2.5 2.5]; % important to keep original EPI vox size
par.jobname      = 'spm_normalize_epi_mask';

% ACTIVATION
img_a = e.getSerie('tedana_ACTIVATION').getVolume('.*ts').removeEmpty;
tmp_exam = [img_a.exam];
warp_field = tmp_exam.getSerie('anat').getVolume('^y');
job_apply_normalize(warp_field, img_a, par);

% RS
par.jobname = 'spm_normalize_epi_rs';
img_rs = e.getSerie('tedana_RS').getVolume('.*ts').removeEmpty;
tmp_exam_rs = [img_rs.exam];
warp_field = tmp_exam_rs.getSerie('anat').getVolume('^y');
job_apply_normalize(warp_field, img_rs, par);

%% Normalize the Tmean used later for PhysIO Noise ROI
par.auto_add_obj = 1;
img = e.gser('run_ACT').removeEmpty.gvol('bet');
tmp_exam = [img.exam];
y = tmp_exam.gser('anat').gvol('^y');

par.jobname = 'spm_normalize_meanepi';
job_apply_normalize(y, img, par);

%% Normalize the Tmean mask used later for PhysIO Noise ROI
par.auto_add_obj = 1;

%% ACT
%e.gser('ACTIVATION').addVolume('mask','mask',1)
img = e.gser('ACTIVATION').getVolume('mask'); %% do same for RS
tmp_exam = [img.exam];
y = tmp_exam.gser('anat').gvol('^y');

par.jobname = 'spm_normalize_maskepi_ACT';
job_apply_normalize(y, img, par);

%% RS
%e.gser('ACTIVATION').addVolume('mask','mask',1)
img = e.gser('run_RS').getVolume('mask'); %% do same for RS
tmp_exam = [img.exam];
y = tmp_exam.gser('anat').gvol('^y');

par.jobname = 'spm_normalize_maskepi_RS';
job_apply_normalize(y, img, par);

%% add outputs to exam object
e.getSerie('tedana').addVolume('^wts','wts',1);
e.getSerie('tedana').addVolume('^wdn','wdn',1);
e.getSerie('run').addVolume('^wbet.*vtde1.nii','wbet',1);
e.getSerie('anat').addVolume('^wp0','wp0',1);
e.getSerie('run').addVolume('^wbet.*mask.nii$','wmask',1)

%% Smooth TEDANA outputs
%% smooth with 5

clear par
par.redo         = 0;
if CLUSTER
    par.run      = 0;
    par.sge      = 1;
else
    par.run      = 1;
    par.sge      = 0;
end
par.auto_add_obj = 0;
par.smooth       = [5 5 5];
par.prefix       = 's5';
par.jobname      = 'spm_smooth_last_5';

img = e.gser('tedana').gvol('^w').removeEmpty;

job_smooth(img, par);

% save smoothed outputs
e.gser('tedana').addVolume('^s5wts_OC.nii','s5wts',1);
e.gser('tedana').addVolume('^s5wdn_ts_OC.nii','s5wdn',1);


%% smooth with 6

clear par
par.redo         = 0;
if CLUSTER
    par.run      = 0;
    par.sge      = 1;
else
    par.run      = 1;
    par.sge      = 0;
end
par.auto_add_obj = 0;
par.smooth       = [6 6 6];
par.prefix       = 's6';
par.jobname      = 'spm_smooth_last_6';

img = e.gser('tedana').gvol('^w').removeEmpty;

job_smooth(img, par);

% save smoothed outputs

e.gser('tedana').addVolume('^s6wts_OC.nii','s6wts',1);
e.gser('tedana').addVolume('^s6wdn_ts_OC.nii','s6wdn',1);

return

%% coregister WM & CSF on functional (warped_mean_mask)
% used for TAPAS: PhysIO

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
else
    par.run = 1;
    par.sge = 0;
end
par.type    = 'estimate_and_write';

ref = e.getSerie('run');
ref = ref(:,1).getVolume('wbet');
e.getSerie('anat').addVolume('^wp2','wp2',1);
e.getSerie('anat').addVolume('^wp3','wp3',1);
src = e.getSerie('anat').getVolume('^wp2');
oth = e.getSerie('anat').getVolume('^wp3');

job_coregister(src, ref, oth, par);

% %% coregister WM & CSF resliced on voxel size 2.5 on functional
% (warped_mean_mask) not needed if CAT12 does its job and creates the rwp2
% & rwp3 files
% % used for TAPAS: PhysIO
% 
% clear par
% if CLUSTER
%     par.run = 0;
%     par.sge = 1;
% else
%     par.run = 1;
%     par.sge = 0;
% end
% par.type    = 'estimate_and_write';
% 
% ref = e.getSerie('run');
% ref = ref(:,1).getVolume('wbet');
% e.getSerie('anat').addVolume('^wp2.*reslice','wp2_reslice',1);
% e.getSerie('anat').addVolume('^wp3.*reslice','wp3_reslice',1);
% src = e.getSerie('anat').getVolume('^wp2_reslice');
% oth = e.getSerie('anat').getVolume('^wp3_reslice');
% 
% job_coregister(src, ref, oth, par);

%% Check & Save
e.explore
cd (main_dir)
save('e','e')
    
%%

