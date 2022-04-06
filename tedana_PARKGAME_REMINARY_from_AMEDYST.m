clear
clc

par.pct = 0;

main_dir = fullfile('/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek','nifti');

e_PARKGAME = exam(main_dir,'PARKGAME');
e_REMINARY = exam(main_dir,'REMINARY_\w{2}_');

%e = e_PARKGAME + e_REMINARY;
e = exam(main_dir,'2018_07_18_PARKGAMEII_001_NB_18_07_2018_V1');

%% Add series

% T1 anat
e.addSerie('t1mpr_S256_0_8iso_p2$','anat_T1',1)
e.getSerie('anat_T1').addVolume('^s.*p2.nii','s',1)

% Resting State
e.addSerie(        'RS$','run_RS'        ,1)

% Activation fMRI
e.addSerie('ACTIVATION$','run_ACTIVATION',1)

e.getSerie('run').addVolume('^f','f',3)

% [ ec, ei ] = e.removeIncomplete();

e.reorderSeries('path');
e.unzipVolume(par);


%% remove bad subjects pas fait

% fprintf('nSubj before = %d \n', numel(e))
% e =
% e.removeTag({'Subject019','Subject038','Subject053','Subject062','Subject066'});
% can be used if we add a loop to verify number of volumes per subject
%fprintf('nSubj after = %d \n', numel(e))


%% check number of volumes ACTIVATION

p = e.getSerie('run_ACTIVATION').getVolume('^f').print;
p = cellstr(p);
res =regexp(p,'f630');
res = cellfun(@isempty,res);
bad_volumes = char(p(res));

if ~isempty(bad_volumes)
    disp(bad_volumes)
    error('wrong number of volumes')
else
    fprintf('all echos seem to have 630 volumes \n')
end

%% check number of volumes RS

p = e.getSerie('run_RS').getVolume('^f').print;
p = cellstr(p);
res =regexp(p,'f300');
res = cellfun(@isempty,res);
bad_volumes = char(p(res));

if ~isempty(bad_volumes)
    disp(bad_volumes)
    error('wrong number of volumes')
else
    fprintf('all echos seem to have 300 volumes \n')
end


%% Sort echos

clear par

par.fname        = 'meinfo'; % name of the .mat file that will be saved
par.sge          = 0;
par.redo         = 0;
par.run          = 1;
par.fake         = 0;

run_list = e.getSerie('run');

[meinfo, jobs] = job_sort_echos( run_list , par );
meinfo.volume = e.getSerie('run').getVolume('^e');
meinfo.anat   = e.getSerie('anat').getVolume('^s');


%% afni_proc

clear par

% TEDANA recomandations :          tshift,        volreg
% MEICA pipeline        : despike, tshift, align, volreg
% robust for TEDANA     : despike, tshift, align, volreg   BUT it "blurs" the data

par.blocks   = {'despike','tshift','volreg'}; % now codded : despike, tshift, align, volreg
par.seperate = 1;                             % each volume is treated seperatly : useful when runs have different orientations
par.execute  = 1;                             % execute afni_proc.py generated tcsh script file immidatly after the generation

par.run      = 0;
par.sge      = 1;

par.jobname  = 'job_afni_proc_multi_echo';
par.subdir   = 'afni';

par.pct      = 0;
par.redo     = 0;
par.fake     = 0;
par.verbose  = 1;

job_afni_proc_multi_echo( meinfo , par )


return

%% pour check results

%e.getSerie('rs').getVolume('vtde1').removeEmpty.getOne
% 1dplot -volreg dfile.r01.1D ==> realignment paramters,
% rp=[Rx,Ry,Rz,Tx,Ty,TZ] (deg° / mm)
% 1dplot motion_Prestim__S38_RS_enorm.1D ==> sum(diff(drp)^2) ~ similira
% computation as FD (Framewise Displacement)


%% get pp from afni

e.getSerie('run_RS').addVolume('^vtde1.nii','vtde1',1)
e.getSerie('run_RS').addVolume('^vtde2.nii','vtde2',1)
e.getSerie('run_RS').addVolume('^vtde3.nii','vtde3',1)


%% temporal mean

clear par
par.run   = 0;
par.fake  = 0;
par.sge   = 1;
par.fsl_output_format = 'NIFTI';

par.jobname = 'fslmean';
vtde_in  = e.getSerie('run_RS').getVolume('^vtde1').removeEmpty.toJob;
vtde_out = addprefixtofilenames(vtde_in,'mean_');
do_fsl_mean(vtde_in,vtde_out,par);

e.getSerie('rs').addVolume('^mean_vtde1.nii','mean_vtde1',1)


%% bet

clear par
par.run   = 0;
par.fake  = 0;
par.sge   = 1;
par.fsl_output_format = 'NIFTI';

run_dir = e.getSerie('rs').removeEmpty.toJob(0);

job = cell(0);
for iJob = 1 : length(run_dir)
    job{iJob,1} = sprintf('export FSLOUTPUTTYPE=NIFTI; cd %s; /network/lustre/iss02/cenir/software/irm/fsl5/bin/bet %s %s -R -m -f 0.3',...
        run_dir{iJob},...
        'mean_vtde1.nii',...
        'bet_mean_vtde1.nii');
end
par.jobname = 'fslbet';
do_cmd_sge(job,par);

e.getSerie('run_RS').addVolume('^bet_mean_vtde1_mask.nii','mask_vtde1',1)
e.getSerie('run_RS').addVolume('^bet_mean_vtde1.nii'     , 'bet_vtde1',1)



%% Tedena


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

% Cluster
par.sge      = 1;               % for ICM cluster, run the jobs in paralle
par.jobname  = 'job_tedana';
par.walltime = '08:00:00';      % HH:MM:SS
par.mem      = 16000;           % MB
par.sge_queu = 'normal,bigmem'; % use both

job_tedana( meinfo, 'vtd', 'tedana_vtd_mle', 'bet_mean_vtde1_mask.nii', par );



return

