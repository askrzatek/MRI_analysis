% preproc_tedana

clear all
clc

%% Initialisation
par.pct = 0;
main_dir = fullfile('/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek','nifti');

e_PARKGAME = exam(main_dir,'PARKGAME');
%e_REMINARY = exam(main_dir,'REMINARY_\w{2}_');

e = e_PARKGAME; % (3:length(e_PARKGAME)); % choose specific
e.addSerie('ACTIVATION$','run_ACTIVATION',1)
e.addSerie(        'RS$','run_RS'        ,1)

e.getSerie('run').addVolume('^f\d{3}','f',3)

e.addSerie('t1mpr_S256_0_8iso_p2$','anat_T1',1)
e.getSerie('anat_T1').addVolume('^s.*p2.nii','s',1)

e.reorderSeries('path');
e.unzipVolume(par);

dir_func = e.getSerie('run') .toJob;
dir_anat = e.getSerie('anat').toJob(0);


%% vol quantity control
% functional ACTI
p = e.getSerie('run_ACTIVATION').getVolume('^f').print;
p = cellstr(p);
res =regexp(p,'f630');
res = cellfun(@isempty,res);
bad_volumes = char(p(res));

if ~isempty(bad_volumes)
    disp(bad_volumes)
    %error('wrong number of volumes')
else
    fprintf('all echos seem to have 630 volumes \n')
end

% functional RS
p = e.getSerie('run_RS').getVolume('^f').print;
p = cellstr(p);
res =regexp(p,'f300');
res = cellfun(@isempty,res);
bad_volumes = char(p(res));

if ~isempty(bad_volumes)
    disp(bad_volumes)
    %error('wrong number of volumes')
else
    fprintf('all echos seem to have 300 volumes \n')
end

%% sort echos

clear par

par.fname        = 'meinfo'; % name of the .mat file that will be saved
par.sge          = 0;
par.redo         = 0;
par.run          = 1;
par.fake         = 0;

run_list = e.getSerie('run');

[meinfo, jobs] = job_sort_echos( run_list , par )

meinfo.volume = e.getSerie('run').getVolume('^e')
meinfo.anat   = e.getSerie('anat').getVolume('^s')

%% AFNI - vtd

clear par
par.blocks   = {'despike','tshift','volreg'}; % now codded : despike, tshift, align, volreg
par.seperate = 1;                             % each volume is treated seperatly : useful when runs have different orientations
par.execute  = 1;                             % execute afni_proc.py generated tcsh script file immidatly after the generation

par.run      = 1;
par.sge      = 0;

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

[ec_vtd, ei_vtd] = e.removeIncomplete;
e = ec_vtd;

%% temporal mean

clear par
par.run   = 1;
par.fake  = 0;
par.sge   = 0;
par.fsl_output_format = 'NIFTI';

par.jobname = 'fslmean';
vtde_in  = e.getSerie('run').getVolume('^vtde1').removeEmpty.toJob;
vtde_out = addprefixtofilenames(vtde_in,'mean_');
do_fsl_mean(vtde_in,vtde_out,par);

e.getSerie('run').addVolume('^mean_vtde1.nii','mean_vtde1',1)


%% bet (brain extraction)

clear par
par.run   = 1;
par.fake  = 0;
par.sge   = 0;
par.fsl_output_format = 'NIFTI';

run_dir = e.getSerie('run').removeEmpty.toJob(0);

job = cell(0);
for iJob = 1 : length(run_dir)
    job{iJob,1} = sprintf('export FSLOUTPUTTYPE=NIFTI; cd %s; /network/lustre/iss01/cenir/software/irm/fsl5/bin/bet %s %s -R -m -f 0.3',...
        run_dir{iJob},...
        'mean_vtde1.nii',...
        'bet_mean_vtde1.nii');
end
par.jobname = 'fslbet';
do_cmd_sge(job,par);

e.getSerie('run').addVolume('^bet_mean_vtde1_mask.nii','mask_vtde1',1)
e.getSerie('run').addVolume('^bet_mean_vtde1.nii'     , 'bet_vtde1',1)


%% Tedana

% tedana.py main arguments
par.tedpca     = 'mle'; % mle, kundu, kundu-stabilize (tedana default = mle)
par.maxrestart = []; % (tedana default = 10)
par.maxit      = []; % (tedana default = 500)
par.png        = 1;  % tedana will make some PNG files of the components Beta map, for quality checks

% tedana.py other arguments
par.cmd_arg = ''; % Allows you to use all addition arguments not scripted in this job_tedana.m file

% matvol classic options
par.pct      = 0; % Parallel Computing Toolbox, will execute in parallel all the subjects
par.redo     = 0; % overwrite previous files
par.fake     = 0; % do everything exept running
par.verbose  = 2; % 0 : print nothing, 1 : print 2 first and 2 last messages, 2 : print all
par.run      = 1;

% Cluster
if par.sge               % for ICM cluster, run the jobs in paralle
    par.jobname  = 'job_tedana';
    par.walltime = '08:00:00';      % HH:MM:SS
    par.mem      = 16000;           % MB
    par.sge_queu = 'normal,bigmem'; % use both
end


job_tedana( meinfo, 'vtd', 'tedana_vtd_mle', 'bet_mean_vtde1_mask.nii', par );


%return %not sure if necessary

%% tedana report % meinfo is in main_dir and tedana outputs in the subdir of run dir (until change - cannot use tedana_report)

%tedana_report (main_dir);

%% Add new preprocessed data to exam object : dn_ts_OC.nii 

e.addSerie('ACTIVATION','tedana_vtd','tedana_ACTIVATION',1);
e.addSerie('RS','tedana_vtd','tedana_RS',1);
e.getSerie('tedana').addVolume('^dn.*.nii$','dn',1);

[ec_dn, ei_dn] = e.removeIncomplete;
e = ec_dn;

% if all in run directory instead of tedana_vtd_mle
%e(1).getSerie('run').addVolume('^dn.*.nii$','dn',1);

%% Run CAT12 segmentation
clear par
%anat segment
fanat = e.getSerie('anat').getVolume('^s').getPath;

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
par.redo    = 0;

j_segment = job_do_tedana_segmentCAT12(fanat,par);


e.getSerie('anat').addVolume('^y'  ,'y' );
e.getSerie('anat').addVolume('^p0' ,'p0' );

%% Coregister func @ anat
% correction of any former transformations of ^dn

    par.ask = 0;
    par.redo = 1;
    par.pct = 0;
    
    if par.redo==1
        
        e.getSerie('tedana').addVolume('^ts','ts',1);
        
        %% for ACTIVATION vols
        origin = e.getSerie('tedana_ACTIVATION').getVolume('^ts').toJob(0);
        influencer = e.getSerie('run_ACTIVATION').getVolume('^bet_vtde1').toJob(0);
        follower = e.getSerie('tedana_ACTIVATION').getVolume('^dn').toJob(0);
                   
        do_fsl_copy_header(influencer, origin,      par);
        do_fsl_copy_header(follower,   influencer,  par);
        

    %% for RS vols
        origin = e.getSerie('tedana_RS').getVolume('^ts').toJob(0);
        influencer = e.getSerie('run_RS').getVolume('^bet_vtde1').toJob(0);
        follower = e.getSerie('tedana_RS').getVolume('^dn').toJob(0);
       
        par.ask = 0;
        par.pct = 0;
        do_fsl_copy_header(influencer, origin,      par);
        do_fsl_copy_header(follower,   influencer,  par);
    
    end

%% Coregistration of both runs

% ACTIVATION
    clear par
    par.type   = 'estimate';
    par.interp = 1;
    par.prefix = 'r';
    par.sge    = 0;
    par.redo   = 0;
    par.run    = 1;
    par.display= 0;

    ref = e.getSerie('anat').getVolume('^p0').toJob(0);
    src = e.getSerie('run_ACTIVATION').getVolume('^bet_vtde1').toJob(0);
    oth = e.getSerie('tedana_ACTIVATION').getVolume('^dn').toJob(0);

        % include the skip option (?)

    job_coregister(char(src),char(ref),char(oth),par)
    
% RS
    src = e.getSerie('run_RS').getVolume('^bet_vtde1').toJob(0);
    oth = e.getSerie('tedana_RS').getVolume('^dn').toJob(0);

    % include the skip option

    job_coregister(char(src),char(ref),char(oth),par)
    
    %% Normalize

    clear par
    par.redo    = 1;
    par.sge     = 0;
    par.run     = 1;
    par.display = 0;
    par.jobname = 'spm_apply_norm';

    warp_field = e.getSerie('anat').getVolume('^y');
    img = e.getSerie('tedana_RS').getVolume('^dn');

    job_apply_normalize(warp_field,img, par)

    e.getSerie('tedana').addVolume('^wdn','wdn',1);
    e.explore

    save('e','e')
    
    %%

