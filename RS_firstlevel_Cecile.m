%% script RS firstlevel : Cecile Gallea

clear
clc

load e

%% Modeling low frequency fluctuations

% -> four binaries Fourier pairs with 90 degree phase lag at the frequency
% of 0.01, 0.02, 0.04, and0.08 Hz (100, 50,25, and 12.5 s),

t = linspace(0,480,300)';
% 0.01 Hz ou 100s
y1 = square(2*pi*0.01*t+(pi/2));
y2 = square(2*pi*0.01*t);
% 0.02 Hz ou 50s
y3 = square(2*pi*0.02*t+(pi/2));
y4 = square(2*pi*0.02*t);
% 0.04 Hz ou 25s
y5 = square(2*pi*0.04*t+(pi/2));
y6 = square(2*pi*0.04*t);
% 0.08 Hz ou 12.5s
y7 = square(2*pi*0.08*t+(pi/2));
y8 = square(2*pi*0.08*t);


%% Prepare dirs

dir_RS   = e.getSerie('RS$').removeEmpty. toJob(0);
dir_RS   = cellstr(dir_RS(1:length(dir_RS)/2));
dirFunc  = e.getSerie('tedana').removeEmpty.toJob(0);

dirStats = e.getSerie('RS$').mkdir('model','model_1'); % basic model ALFF
dirStats = dirStats(1:length(dirStats)/2);
%% Make symbolic links from tedana_vtd_mle dir to run dir based on job_meica_afni symbolic link creation 
addpath /home/anna.skrzatek/MRI_analysis/
% %par.subdir        = 'tedana_vtd_mle';
% par.subdir        = 'tedana009a1_vtd';
% par.sge = 0;
% par.run = 1;
% 
% par.warp_file_reg = '^s6wts';
% job_symbolic_child_to_parent(dir_RS, par);
% 
%% make a symbolic link of rp_spm.txt and of multiple_regressors to dir_func
clear par
par.subdir = 'wts';
par.regfile_regex = 'multiple_regressors.txt';
for i= 1:length(e)
    wd = e(i).getSerie('RS').path;
    reg_dir = char(get_subdir_regex(wd, par.subdir));
    A_src = fullfile(reg_dir, par.regfile_regex);
    regfile_out = sprintf('%s_%s',par.subdir,par.regfile_regex);
    A_dst = fullfile(wd, regfile_out);
    
    par.redo = 0;
    par.verbose = 2;
    %par.run = 0;
    par.run = 1;
    par.jobname = sprintf('job_symbolic_link');
    %par.jobname = sprintf('%s_%s_%s', 'job_symbolic_link', wd(end-46:end-16), par.subdir);
    [job_session(i)] = r_movefile(A_src, A_dst, 'linkn', par);
    job = [job_session];
    
end
%job = do_cmd_sge(job, par);

%%

% e.getSerie('RS').mkdir('model','model_CBVI_L');
% e.getSerie('RS').mkdir('model','model_CBVI_R');
% e.getSerie('RS').mkdir('model','model_CBVIII_L');
% e.getSerie('RS').mkdir('model','model_CBVIII_R');


% e.addSerie('Postsham', 'S\d\d_RS$', 'model', 'model_1', 'm1_Postsham', 1 );
% e.addSerie('Poststim', 'S\d\d_RS$', 'model', 'model_1', 'm1_Poststim', 1 );
% e.addSerie('Presham', 'S\d\d_RS$', 'model', 'model_1', 'm1_Presham', 1 );
% e.addSerie('Prestim', 'S\d\d_RS$','model', 'model_1', 'm1_Prestim', 1 );
% dirStats = e.getSerie('m1').toJob;
%dirStats = e.mkdir('model','model_1');


for iSubj = 1 : length(dirFunc)
    for iRun = 1
    %run_dir = dirFunc{iSubj}{iRun};
   
    u1{iSubj}{iRun} = struct('name', 'LFF_0.01Hz_1', 'val',y1);
    user_reg1 = u1';
    u2{iSubj}{iRun} = struct('name', 'LFF_0.01Hz_2', 'val',y2);
    user_reg2 = u2';
    u3{iSubj}{iRun} = struct('name', 'LFF_0.02Hz_1', 'val',y3);
    user_reg3 = u3';
    u4{iSubj}{iRun} = struct('name', 'LFF_0.02Hz_2', 'val',y4);
    user_reg4 = u4';
    u5{iSubj}{iRun} = struct('name', 'LFF_0.04Hz_1', 'val',y5);
    user_reg5 = u5';
    u6{iSubj}{iRun} = struct('name', 'LFF_0.04Hz_2', 'val',y6);
    user_reg6 = u6';
    u7{iSubj}{iRun} = struct('name', 'LFF_0.08Hz_1', 'val',y7);
    user_reg7 = u7';
    u8{iSubj}{iRun} = struct('name', 'LFF_0.08Hz_2', 'val',y8);
    user_reg8 = u8';
   
     end % iRun
end % iSubj

clear y1 y2 y3 y4 y5 y6 y7 y8 t Matrix Mat_sub u1 u2 u3 u4 u5 u6 u7 u8


%% Parameters

par.TR = 1.6;

par.file_reg = '^s6wts.*nii';
par.rp       = 1;
par.rp_regex = 'wts_multiple.*txt';

% Masking
par.mask_thr = 0.8; % spm default option
par.mask     =  {}; % cell(char) of the path for the mask of EACH model : N models means N paths
par.cvi      = 'AR(1)'; % 'AR(1)' / 'FAST' / 'none'

% Regressors
%-----------
% multilevel_cells(struct) for used defined regressors : they will NOT be convolved
par.user_regressor = {};
% multilevel_cells(char  ) for used defined regressors : they will NOT be convolved
% The regressors in the file will be concatenated with rp_*.txt
par.file_regressor = addsuffixtofilenames(dir_RS, '/wts_multiple_regressors.txt');

par.jobname  = 'spm_glm';
par.walltime = '04:00:00';

par.sge   = 0;
par.run      = 0;
par.display  = 1;
par.redo     = 0;


% %% tests
% dirFunc = dirFunc(1);
% dirFunc = dirFunc(1);
% dirStats = dirStats(1);


%% Job
%jobs = job_first_level_specify(dirFunc,dirStats, par)


nrSubject=length(dirFunc);


skip = [];

for subj = 1:nrSubject


    spm_file = char(addsuffixtofilenames(dirStats(subj),'SPM.mat'));
    if ~par.redo   &&  exist(spm_file,'file')
        skip = [skip subj];
        fprintf('[%s]: skiping subj %d because %s exist \n',mfilename,subj,spm_file);
    else


        subjectRuns = get_subdir_regex_files(dirFunc{subj},par.file_reg);
        unzip_volume(subjectRuns);
        subjectRuns = get_subdir_regex_files(dirFunc{subj},par.file_reg,struct('verbose',0));
        if par.rp
            fileRP = get_subdir_regex_files(dir_RS{subj},par.rp_regex);
        end


        for run = 1:length(subjectRuns)
            currentRun = cellstr(subjectRuns{run}) ;
            clear allVolumes

            if length(currentRun) == 1 %4D file
                allVolumes = spm_select('expand',currentRun);
            else
                allVolumes = currentRun;
            end
            jobs{subj}.spm.stats.fmri_spec.sess(run).scans = allVolumes; %#ok<*AGROW>
            jobs{subj}.spm.stats.fmri_spec.sess(run).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            jobs{subj}.spm.stats.fmri_spec.sess(run).multi = {''};


            if par.rp && isempty(par.file_regressor)
                jobs{subj}.spm.stats.fmri_spec.sess(run).multi_reg = fileRP(run);
            elseif ~par.rp && ~isempty(par.file_regressor)
                jobs{subj}.spm.stats.fmri_spec.sess(run).multi_reg = par.file_regressor(subj);
            elseif par.rp && ~isempty(par.file_regressor)
                jobs{subj}.spm.stats.fmri_spec.sess(run).multi_reg = [par.file_regressor{subj}(run) ; fileRP(run)];
            else
                jobs{subj}.spm.stats.fmri_spec.sess(run).multi_reg = {''};
            end

            jobs{subj}.spm.stats.fmri_spec.sess(run).regress = struct('name', {}, 'val', {});
            jobs{subj}.spm.stats.fmri_spec.sess(run).regress(1) = user_reg1{subj}{run};
            jobs{subj}.spm.stats.fmri_spec.sess(run).regress(2) = user_reg2{subj}{run};
            jobs{subj}.spm.stats.fmri_spec.sess(run).regress(3) = user_reg3{subj}{run};
            jobs{subj}.spm.stats.fmri_spec.sess(run).regress(4) = user_reg4{subj}{run};
            jobs{subj}.spm.stats.fmri_spec.sess(run).regress(5) = user_reg5{subj}{run};
            jobs{subj}.spm.stats.fmri_spec.sess(run).regress(6) = user_reg6{subj}{run};
            jobs{subj}.spm.stats.fmri_spec.sess(run).regress(7) = user_reg7{subj}{run};
            jobs{subj}.spm.stats.fmri_spec.sess(run).regress(8) = user_reg8{subj}{run};


            if ~isempty(par.user_regressor)
                jobs{subj}.spm.stats.fmri_spec.sess(run).regress = par.user_regressor{subj}{run};
            end
            jobs{subj}.spm.stats.fmri_spec.sess(run).hpf = 128;

        end % run

        jobs{subj}.spm.stats.fmri_spec.timing.units = 'secs';
        jobs{subj}.spm.stats.fmri_spec.timing.RT = par.TR;
        jobs{subj}.spm.stats.fmri_spec.timing.fmri_t = 16;
        jobs{subj}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

        jobs{subj}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        jobs{subj}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        jobs{subj}.spm.stats.fmri_spec.volt = 1;
        jobs{subj}.spm.stats.fmri_spec.global = 'None';
        jobs{subj}.spm.stats.fmri_spec.mthresh = par.mask_thr;
        if isempty(par.mask)
            jobs{subj}.spm.stats.fmri_spec.mask = {''};
        else
            jobs{subj}.spm.stats.fmri_spec.mask = par.mask(subj);
        end
        jobs{subj}.spm.stats.fmri_spec.cvi = par.cvi;

    end % SPM.mat exists ?

    jobs{subj}.spm.stats.fmri_spec.dir = dirStats(subj);

end


%% Other routines
% verifs
main_dir = '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test';
% subj_dir = gdir(main_dir,'^Subj|^___S')
subj_dir = gdir(main_dir,'PARKGAME.*[a,c]$')
SessDir = gdir(subj_dir,'.*RS$');
StatDir = gdir(SessDir,'model','model_1');
regressors_test = gfile(SessDir,'wts_multiple_regressors');

par.run=0;
par.sge =1;
par.display = 0;
par.jobname = 'spm_first_level_spec_RS'

[ jobs ] = job_ending_rountines( jobs, skip, par );
%% Estimate

% main_dir = '/network/lustre/iss01/cenir/analyse/irm/studies/AMEDYST/RS/nifti';
% subj_dir = gdir(main_dir,'^Subj|^')
%
% spm_file = gfile(dirStats(subj),'SPM.mat'))
%
% Modele_Sess_dir = gdir(subj_dir,'^P','RS$','^modele','model_1');
% SPM_Sess_file = gfile(Modele_Sess_dir,'SPM.mat');

%fspm = addsuffixtofilenames( dirStats, 'SPM.mat');
fspm = addsuffixtofilenames(StatDir, 'SPM.mat');

clear par
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.jobname  = 'spm_first_level_est_RS';
job_first_level_estimate(fspm,par)

%% Contrast estimation

clear par
par.sge=1;
par.sge_queu = 'normal,bigmem';
par.jobname = 'spm_first_level_con_RS';


EffectsOfInterest= [1 0 0 0 0 0 0 0; ...
                    0 1 0 0 0 0 0 0; ...
                    0 0 1 0 0 0 0 0; ...
                    0 0 0 1 0 0 0 0; ...
                    0 0 0 0 1 0 0 0; ...
                    0 0 0 0 0 1 0 0; ...
                    0 0 0 0 0 0 1 0; ...
                    0 0 0 0 0 0 0 1];
PositiveEffectALFF = [1 1 1 1 1 1 1 1];
NegativeEffectALFF = [1 1 1 1 1 1 1 1];
PositiveLowFreq   = [1 1 1 1 0 0 0 0];
PositiveHighFreq   = [0 0 0 0 1 1 1 1];

contrast.names={'Effects Of Interest','PositiveEffects','NegativeEffects','PositiveLowFreq','PositiveHighFreq'};
contrast.values={EffectsOfInterest,PositiveEffectALFF,NegativeEffectALFF,PositiveLowFreq,PositiveHighFreq};
contrast.types={'F','T','T','T','T'};

par.delete_previous=1

j=job_first_level_contrast(fspm,contrast,par)
