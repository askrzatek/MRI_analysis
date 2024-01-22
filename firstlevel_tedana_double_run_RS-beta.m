%% adapted from script RS firstlevel : Cecile Gallea 
%% Init

clear
clc

addpath /home/anna.skrzatek/matvol/SPM/firstlevel/
addpath('/home/anna.skrzatek/MRI_analysis/')

main_dir = '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test';
cd(main_dir)

%addpath '/network/lustre/iss02/cenir/analyse/irm/users/cecile.gallea/ASYA/asyasuit/'

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

model_name = {'full_RS_sts_tapas_doublerun_resliced'};%, 'smodel_ts_tapas', 'smodel_dn_tapas'};
mkdir(main_dir,model_name{1});
double_model_dir = fullfile(main_dir,model_name{1});

%% fetch Input dirs
RSinput_dir = fullfile(main_dir,'/firstlevel_RS')
cd (RSinput_dir)

patient_regex = {'PARKGAMEII.*NB.*_a','PARKGAMEII.*BM.*_a','PARKGAMEII.*SM.*_c','PARKGAMEII.*SD.*_a','PARKGAMEII.*JR.*_a','PARKGAMEII.*LJ.*_c','PARKGAMEII.*CA.*_a','PARKGAMEII.*PC.*_c','PARKGAMEII.*DD.*','PARKGAMEII.*KM.*_a','PARKGAMEII.*PD.*_a','PARKGAMEII.*CK.*_c','PARKGAMEII.*BF.*_c','PARKGAMEII.*SB.*_a','PARKGAMEII_052.*HJ.*_c'}; %,'PARKGAMEII.*LM.*_c'};
%patient_regex = {'PARKGAMEII_001_NB_a','PARKGAMEII_002_BM_a','PARKGAMEII_003_SM_c','PARKGAMEII_007_SD_a','PARKGAMEII_008_JR_a','PARKGAMEII_023_LJ_c','PARKGAMEII_025_CA_a','PARKGAMEII_028_PC_c','PARKGAMEII_033_DD','PARKGAMEII_039_KM_a','PARKGAMEII_043_PD_a','PARKGAMEII_044_CK_c','PARKGAMEII_047_BF_c','PARKGAMEII_048_SB_a','PARKGAMEII_052_HJ_c'}; 
% patient_regex = {'PARKGAMEII.*009_HJ','PARKGAMEII.*013_RP','PARKGAMEII.*027_OR','PARKGAMEII.*046_HJ','PARKGAMEII.*053_LM}; %exclu

ROIs = {'Caudate_L','Caudate_R','Cereb3_L','Cereb3_R','Cereb4_5_L','Cereb4_5_R','Cereb6_L','Cereb6_R','Cereb7b_L','Cereb7b_R','Cereb8_L','Cereb8_R','Cereb9_L','Cereb9_R','Cereb10_L','Cereb10_R','Cereb_Crus1_L','Cereb_Crus1_R','Cereb_Crus2_L','Cereb_Crus2_R','Cuneus_L','Cuneus_R','Insula_L','Insula_R','Pallidum_L','Pallidum_R','Parietal_Inf_L','Parietal_Inf_R','Parietal_Sup_L','Parietal_Sup_R','PCC','Postcentral_L','Postcentral_R','PPN_L','PPN_R','Precentral_L','Precentral_R','Precuneus_L','Precuneus_R','Putamen_L','Putamen_R','SMA_Ant_L','SMA_Ant_R','SMA_Post_L','SMA_Post_R','Thalamus_L','Thalamus_R','Vermis1_2','Vermis3','Vermis4_5','Vermis6','Vermis7','Vermis8','Vermis9','Vermis10'};
%    ROIs = {'SMA_face_L','SMA_face_R','SMA_foot_L','SMA_foot_R','SMA_hand_L','SMA_hand_R'};

% patient_regex = {'PARKGAMEII.*NB.*_a'}
% ROIs = {'Caudate_L'};

dirFonc = cell(length(patient_regex),length(ROIs));
RS_V1 = exam(main_dir, 'RR');
RS_V2 = RS_V1;

i = 1;
for ip = 1 : length(patient_regex)    
    clear esuj
    esuj = exam(RSinput_dir,patient_regex{ip});
    Stat_dir = gpath(esuj);
    
    e = exam(main_dir,patient_regex{ip});
    e.addSerie('RS$','run_RS',1);
    reg_src_dir1 = e(1).getSerie('run_RS') .toJob;
    reg_src_dir2 = e(2).getSerie('run_RS') .toJob;
        
    for ir = 1 : length(ROIs)
        esuj.addSerie(sprintf('%s',ROIs{ir}),sprintf('%s',ROIs{ir}),2);
        if length(esuj.getSerie(ROIs{ir})) == 2
            dirFonc(ip,ir) = esuj.getSerie(sprintf('%s',ROIs{ir})) .toJob;
        end


% symbolic link multiple regressors in input dirs (firstlevel/subject_name/ROI_name)
        par.subdir = 'wts';
        par.regfile_regex = 'multiple_regressors.txt';
        regfile_out1 = sprintf('%s_%s_V1.txt',par.subdir,par.regfile_regex(1:end-4));
        regfile_out2 = sprintf('%s_%s_V2.txt',par.subdir,par.regfile_regex(1:end-4));

        A_src1 = fullfile(reg_src_dir1{:}, sprintf('%s_%s',par.subdir,par.regfile_regex));
        A_src2 = fullfile(reg_src_dir2{:}, sprintf('%s_%s',par.subdir,par.regfile_regex));
        
        %% needs to be changed to the generic subject directory
        A_dst1 = fullfile(Stat_dir, regfile_out1);
        A_dst2 = fullfile(Stat_dir, regfile_out2);

        par.redo = 0;
        par.verbose = 2;
        par.run = 0;
        par.run = 1;
        par.jobname = sprintf('job_symbolic_link');
        %par.jobname = sprintf('%s_%s_%s', 'job_symbolic_link', wd(end-46:end-16), par.subdir);
        [job_session(i)] = r_movefile(A_src1, A_dst1, 'linkn', par);
        job = [job_session];

        [job_session(i)] = r_movefile(A_src2, A_dst2, 'linkn', par);
        job = [job_session];

% symbolic link wbet-mask in input dirs (firstlevel/subject_name/ROI_name)
        clear par
        par.regfile_regex = 'wbet_Tmean_vtde1_mask.nii';
        regfile_out1 = sprintf('%s_V1.nii',par.regfile_regex(1:end-4));
        regfile_out2 = sprintf('%s_V2.nii',par.regfile_regex(1:end-4));

        A_src1 = fullfile(reg_src_dir1{:}, par.regfile_regex);
        A_src2 = fullfile(reg_src_dir2{:}, par.regfile_regex);
        A_dst1 = fullfile(Stat_dir, regfile_out1);
        A_dst2 = fullfile(Stat_dir, regfile_out2);

            
        par.redo = 0;
        par.verbose = 2;
        par.run = 1;
        par.jobname = sprintf('job_symbolic_link');
        [job_session(i)] = r_movefile(A_src1, A_dst1, 'linkn', par);
        job = [job_session];

        [job_session(i)] = r_movefile(A_src2, A_dst2, 'linkn', par);
        job = [job_session];

% symbolic link s6wts in input dirs (firstlevel/subject_name/ROI_name)
        clear par
        par.regfile_regex = '^s6wts'
        A_src1 = get_subdir_regex_files(reg_src_dir1{:}, par.regfile_regex);
        regfile_out1 = char(A_src1);
        regfile_out1 = sprintf('%s_V1.nii',regfile_out1(end-11:end-4));
        A_src2 = get_subdir_regex_files(reg_src_dir2{:}, par.regfile_regex);
        regfile_out2 = char(A_src1);
        regfile_out2 = sprintf('%s_V2.nii',regfile_out2(end-11:end-4));
        A_dst1 = fullfile(Stat_dir, regfile_out1);
        A_dst2 = fullfile(Stat_dir, regfile_out2);

        par.redo = 0;
        par.verbose = 2;
        par.run = 1;
        par.jobname = sprintf('job_symbolic_link');
        [job_session(i)] = r_movefile(A_src1, A_dst1, 'linkn', par);
        job = [job_session];

        [job_session(i)] = r_movefile(A_src2, A_dst2, 'linkn', par);
        job = [job_session];
        
    end
    
    for iRun = 1 : length(e)
        %run_dir = dirFunc{iSubj}{iRun};
       
        u1{ip}{iRun} = struct('name', 'LFF_0.01Hz_1', 'val',y1);
        user_reg1 = u1';
        u2{ip}{iRun} = struct('name', 'LFF_0.01Hz_2', 'val',y2);
        user_reg2 = u2';
        u3{ip}{iRun} = struct('name', 'LFF_0.02Hz_1', 'val',y3);
        user_reg3 = u3';
        u4{ip}{iRun} = struct('name', 'LFF_0.02Hz_2', 'val',y4);
        user_reg4 = u4';
        u5{ip}{iRun} = struct('name', 'LFF_0.04Hz_1', 'val',y5);
        user_reg5 = u5';
        u6{ip}{iRun} = struct('name', 'LFF_0.04Hz_2', 'val',y6);
        user_reg6 = u6';
        u7{ip}{iRun} = struct('name', 'LFF_0.08Hz_1', 'val',y7);
        user_reg7 = u7';
        u8{ip}{iRun} = struct('name', 'LFF_0.08Hz_2', 'val',y8);
        user_reg8 = u8';
   
     end % iRun

    e.getSerie('RS').addVolume('^s6wts','s6wts',1);
    e.getSerie('RS').addVolume('wts_multiple_regressors','wts_rp',1);
    e.getSerie('RS').addVolume('^wbet.*mask','wbet_mask',1);
    
    Stat_all{ip} = Stat_dir;

    RS_V1= RS_V1 + e(1);
    RS_V2= RS_V2 + e(2);

end

clear y1 y2 y3 y4 y5 y6 y7 y8 t Matrix Mat_sub u1 u2 u3 u4 u5 u6 u7 u8
%clear reg_

RS_V1 = RS_V1(2:end);
RS_V2 = RS_V2(2:end);
RS_all = RS_V1 + RS_V2;
RS_dir = RS_all.getSerie('RS').removeEmpty .toJob(0)


%% Parameters
clear par

par.TR = 1.6;

par.file_reg = '^s6wts.*nii';
par.rp       = 1;
par.rp_regex = 'wts_multiple.*txt';

% Masking
par.mask_thr = 0.1; % spm default option
%par.mask     =  {}; % cell(char) of the path for the mask of EACH model : N models means N paths
par.mask     =  get_subdir_regex_files(Stat_all,'wbet.*mask_V1'); % cell(char) of the path for the mask of EACH model : N models means N paths
par.cvi      = 'AR(1)'; % 'AR(1)' / 'FAST' / 'none'

% Regressors
%-----------
% multilevel_cells(struct) for used defined regressors : they will NOT be convolved
par.user_regressor = {};
% multilevel_cells(char  ) for used defined regressors : they will NOT be convolved
% The regressors in the file will be concatenated with rp_*.txt
par.file_regressor = get_subdir_regex_files(Stat_all,'wts_multiple');

par.jobname  = 'spm_glm_rs_wbet';
par.walltime = '04:00:00';

par.sge   = 0;
par.run      = 0;
par.display  = 1;
par.redo     = 0;


%% JOB Cecile

nrSubject=length(RS_dir)/2;

skip = [];

for subj = 1:length(Stat_all)

    spm_file = char(addsuffixtofilenames(Stat_all{subj},'SPM.mat'));
    if ~par.redo   &&  exist(spm_file,'file')
        skip = [skip subj];
        fprintf('[%s]: skiping subj %d because %s exist \n',mfilename,subj,spm_file);
    else

    clear subjectRuns
    
        subjectRuns = gfile(char(Stat_all{subj}{:}),par.file_reg);
        subjectRuns = cellstr(subjectRuns{:});
        unzip_volume(subjectRuns);
        
%        subjectRuns = get_subdir_regex_files(Stat_all{subj}{:},par.file_reg,struct('verbose',0));
%        subjectRuns = cellstr(subjectRuns{:});
        if par.rp
            fileRP = gfile(Stat_all{subj}{:},par.rp_regex);
            fileRP = cellstr(fileRP{:});
        end

        for run = 1:length(subjectRuns)
            currentRun = cellstr(subjectRuns{run}) ;
            currentRP = fileRP(run);
            
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
                jobs{subj}.spm.stats.fmri_spec.sess(run).multi_reg = currentRP;
            elseif ~par.rp && ~isempty(par.file_regressor)
                jobs{subj}.spm.stats.fmri_spec.sess(run).multi_reg = currentRP;
            elseif par.rp && ~isempty(par.file_regressor)
                jobs{subj}.spm.stats.fmri_spec.sess(run).multi_reg = currentRP;
%                 jobs{subj}.spm.stats.fmri_spec.sess(run).multi_reg = [par.file_regressor{subj}(run) ; fileRP(run)];
%                 jobs{subj}.spm.stats.fmri_spec.sess(run).multi_reg = [par.file_regressor{subj} ; currentRP];
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

    jobs{subj}.spm.stats.fmri_spec.dir = Stat_all{subj};

end

%% Other routines
%subj_dir = gdir(main_dir,'PARKGAME.*[a,c]$')

par.run = 0;
par.sge = 1;
par.display = 0;
%par.jobname = 'spm_first_level_spec_RS_wbet';

[ jobs ] = job_ending_rountines( jobs, skip, par );

%% Estimate model
fspm = addsuffixtofilenames(Stat_all, 'SPM.mat');

clear par
%par.run = 1;
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.jobname  = 'spm_first_level_est_RS_wbet';
job_first_level_estimate(fspm,par)

%% Contrast : definition

%% Getting the number of regressors per session - size of the first session contrast //nb of zeros to add

regcountV1 = [];
regcountV2 = [];
regressors = gfile(Stat_all,'wts_multiple_regressors');
for subj = 1 : length(patient_regex)
    regpaths = cellstr(regressors{subj});
    regfileV1 = load(regpaths{1});
    [nvols1 nreg1] = size(regfileV1);
    regcountV1(subj)= nreg1;
    regconV1{subj} = zeros(8,nreg1);
    regfileV2 = load(regpaths{2});
    [nvols2 nreg2] = size(regfileV2);
    regcountV2(subj)= nreg2;
    regconV2{subj} = zeros(8,nreg2);
end
%%

clear par
par.sge = 0;
par.run = 1;
par.sge_queu = 'normal,bigmem';
par.jobname = 'spm_first_level_dble_con_RS_wbet';

par.sessrep = 'none';

if par.sessrep == 'none'
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

    
    for subj = 1: length(fspm)
        
        %% F-con
        EffectsOfInterestV1    = horzcat(EffectsOfInterest, regconV1{subj},zeros(8,8), regconV2{subj});
        EffectsOfInterestV2    = horzcat(zeros(8,8), regconV1{subj}, EffectsOfInterest, regconV2{subj});
        EffectsOfInterest_all  = vertcat(EffectsOfInterestV1, EffectsOfInterestV2);
        EffectsOfInterestV1_V2 = vertcat(EffectsOfInterestV1, -EffectsOfInterestV2);
        EffectsOfInterestV2_V1 = vertcat(-EffectsOfInterestV1, EffectsOfInterestV2);
        
        %% t-con
        Positive_all            = horzcat(PositiveEffectALFF,regconV1{subj}(1,:),PositiveEffectALFF,regconV2{subj}(1,:));
        Negative_all            = horzcat(NegativeEffectALFF,regconV1{subj}(1,:),NegativeEffectALFF,regconV2{subj}(1,:));
        PositiveLowFreq_all     = horzcat(PositiveLowFreq,regconV1{subj}(1,:),PositiveLowFreq,regconV2{subj}(1,:));
        PositiveHighFreq_all    = horzcat(PositiveHighFreq,regconV1{subj}(1,:),PositiveHighFreq,regconV2{subj}(1,:));        
        
        PositiveV1          = horzcat(PositiveEffectALFF,regconV1{subj}(1,:),-PositiveEffectALFF,regconV2{subj}(1,:));
        NegativeV1          = horzcat(NegativeEffectALFF,regconV1{subj}(1,:),-NegativeEffectALFF,regconV2{subj}(1,:));
        PositiveLowFreqV1   = horzcat(PositiveLowFreq,regconV1{subj}(1,:),-PositiveLowFreq,regconV2{subj}(1,:));
        PositiveHighFreqV1  = horzcat(PositiveHighFreq,regconV1{subj}(1,:),PositiveHighFreq,regconV2{subj}(1,:));        
        
        PositiveV2          = horzcat(-PositiveEffectALFF,regconV1{subj}(1,:),PositiveEffectALFF,regconV2{subj}(1,:));
        NegativeV2          = horzcat(-NegativeEffectALFF,regconV1{subj}(1,:),NegativeEffectALFF,regconV2{subj}(1,:));
        PositiveLowFreqV2   = horzcat(-PositiveLowFreq,regconV1{subj}(1,:),PositiveLowFreq,regconV2{subj}(1,:));
        PositiveHighFreqV2  = horzcat(-PositiveHighFreq,regconV1{subj}(1,:),PositiveHighFreq,regconV2{subj}(1,:));       
        
        %%
    



        contrast.names={'Effects Of Interest V1','Effects Of Interest V2','Effects Of Interest all','Effects Of Interest V1>V2','Effects Of Interest V1<V2','PositiveEffects all','NegativeEffects all','PositiveLowFreq all','PositiveHighFreq all','PositiveEffects V1>V2','NegativeEffects V1>V2','PositiveLowFreq V1>V2','PositiveHighFreq V1>V2','PositiveEffects V1<V2','NegativeEffects V1<V2','PositiveLowFreq V1<V2','PositiveHighFreq V1<V2'};
        contrast.values={EffectsOfInterestV1,EffectsOfInterestV2,EffectsOfInterest_all,EffectsOfInterestV1_V2,EffectsOfInterestV2_V1, Positive_all,Negative_all,PositiveLowFreq_all,PositiveHighFreq_all, PositiveV1,NegativeV1,PositiveLowFreqV1,PositiveHighFreqV1, PositiveV2,NegativeV2,PositiveLowFreqV2,PositiveHighFreqV2};
        contrast.types={'F','F','F','F','F','T','T','T','T','T','T','T','T','T','T','T','T'};

        par.delete_previous=1

        j=job_first_level_contrast(fspm(subj),contrast,par)
    end
end

%% VOI extraction
% fspm %check
% fmask = gpath(e.gser('run_RS').gvol('wmask'));

clear par
par.roi_dir = sprintf('%s/ROI_pariet_mot_premot_cereb_BG_PPN',main_dir);
fileROI = cellstr(char(gfile(par.roi_dir,'.*')));

for subj = 1 : length(fspm)
    fmask = gfile(char(Stat_all{subj}{:}),'^mask.nii');
    for iRun = 1 : 2
        par.run = iRun
        
        par.jobname = 'spm_voi_ts_extract_atlas_wbet';
        spm_job_voi(fspm(subj), fmask, par);

        par.jobname = 'spm_voi_PCC_wbet';
        spm_job_voi_PCC(fspm(subj),fmask,par);
        
    end
end

pccmodeldir = get_subdir_regex(RSinput_dir,'PARKGAME.*[a,c]$');
filePCC = cellstr(char(gfile(pccmodeldir,'.*PCC.*.mat')));

char(fileROI)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create output directories architecture
patient_list = {'PARKGAMEII_001_NB_a','PARKGAMEII_002_BM_a','PARKGAMEII_003_SM_c','PARKGAMEII_007_SD_a','PARKGAMEII_008_JR_a','PARKGAMEII_023_LJ_c','PARKGAMEII_025_CA_a','PARKGAMEII_028_PC_c','PARKGAMEII_033_DD','PARKGAMEII_039_KM_a','PARKGAMEII_043_PD_a','PARKGAMEII_044_CK_c','PARKGAMEII_047_BF_c','PARKGAMEII_048_SB_a','PARKGAMEII_052_HJ_c'}; 

%%

contrast_PPI = {
    'Effect of Interest'   1
    };

% model_dir = gdir(subj_dir,'^P','RS$','LFF_BOX_glm');
%models    = fspm;

nRun   = 2;
nSubj = length(Stat_all);
nROI  = length(fileROI)+1;

%?
nCon  = size(contrast_PPI,1);
%?

volume_file = cell(length(Stat_all),2);
for subj = 1 : length(Stat_all)
    vol = gfile(Stat_all{subj},'^s6wts');
    vol = cellstr(vol{:});
    volume_file{subj,1} = vol(1);
    volume_file{subj,2} = vol(2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verif VOI
%% verification creation VOIs

voi_file = cell(nSubj,nROI,nRun);
for subj = 1 : nSubj
    mkdir(double_model_dir, patient_list{subj});
    patients_dir{subj} = get_subdir_regex(double_model_dir, patient_list{subj});
    for ir = 1 : nROI-1
       
       [~,roi_name] = fileparts(fileROI{ir});
       roi_name = roi_name(1:end-10);
       
       mkdir(char(patients_dir{subj}),roi_name);
       outdirs{subj,ir} =  get_subdir_regex(patients_dir{subj},roi_name);%% outdirs
       
       for iRun = 1 : nRun
           voi_file{subj,ir,iRun} = gfile(Stat_all{subj},sprintf('VOI_%s.*_%d.mat',roi_name,iRun));
       end
    end
    [~,roi_name] = fileparts(filePCC{subj});
    roi_name = roi_name(5:end-2);
    
    mkdir(char(patients_dir{subj}),roi_name);
    outdirs{subj,nROI} =  get_subdir_regex(patients_dir{subj},roi_name);%% outdirs
       
    for iRun = 1 : nRun
       voi_file{subj,nROI,iRun} = gfile(Stat_all{subj},sprintf('VOI_%s.*_%d.mat',roi_name,iRun));
    end   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare jobs

% Onsets identical for everybody (resting state)
TR    = 1.6;
nbvol = 300;
jobs = cell(nSubj,nROI);


for subj = 1 : nSubj
    stat_dir = Stat_all{subj};
           
    for ir = 1 : nROI
        if ir == nROI
            roi_name = 'PCC';
        else
            [~,roi_name] = fileparts(fileROI{ir}); % not aplied since we create one VOI individually from a sphere at coordinates (PCC)
        end
       
        for iRun = 1 : nRun
           
            jobs{subj,ir}.spm.stats.fmri_spec.dir = outdirs{subj,ir}; %'PCC')); %VOI__rCerebellum_6_Left__run1_1
            jobs{subj,ir}.spm.stats.fmri_spec.timing.units = 'secs';
            jobs{subj,ir}.spm.stats.fmri_spec.timing.RT = 1.6;
            jobs{subj,ir}.spm.stats.fmri_spec.timing.fmri_t = 16;
            jobs{subj,ir}.spm.stats.fmri_spec.sess(iRun).timing.fmri_t0 = 8;
           
            % VOI
            voi_file = fullfile(Stat_all{subj},sprintf('VOI_%s_%d.mat', roi_name,iRun)); %'PCC')); %
            voi = load(char(voi_file));
            jobs{subj,ir}.spm.stats.fmri_spec.sess(iRun).regress.name = 'Y';
            jobs{subj,ir}.spm.stats.fmri_spec.sess(iRun).regress.val = voi.Y;
           
            % Volumes
            jobs{subj,ir}.spm.stats.fmri_spec.sess(iRun).scans = spm_select('expand',cellstr(volume_file{subj,iRun}));

            % Other
            jobs{subj,ir}.spm.stats.fmri_spec.sess(iRun).multi = {''}; % we don't use .mat file for the onsets
            jobs{subj,ir}.spm.stats.fmri_spec.sess(iRun).hpf = 128;
           
        end % iRun
       
        jobs{subj,ir}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        jobs{subj,ir}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        jobs{subj,ir}.spm.stats.fmri_spec.volt = 1;
        jobs{subj,ir}.spm.stats.fmri_spec.global = 'None';
%         jobs{subj,ir,iRun}.spm.stats.fmri_spec.mthresh = 0; % because already skullstriped before TEDANA
%         jobs{subj,ir,iRun}.spm.stats.fmri_spec.mask = {''};
        jobs{subj,ir}.spm.stats.fmri_spec.mthresh = 0.1; % consistent with previous mask parameters -> according to Benoit
        jobs{subj,ir}.spm.stats.fmri_spec.mask = fullfile(Stat_all{subj},'mask.nii'); % masks from the firstlevel of model_1 or model_2 

        jobs{subj,ir}.spm.stats.fmri_spec.cvi = 'AR(1)';
       
       
       
    end % ir
   
end % subj

par.sge = 0;
par.run = 1;
par.redo = 1;
par.jobname = 'RS_voi_doublerun_model_spec';

job_ending_rountines(jobs,[],par);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  Estimate
cd(double_model_dir)
clear subj ir
for subj = 1 : nSubj
    for ir = 1: nROI
        fspm = gfile(outdirs{subj,ir},'SPM.mat');
    
        clear par
        par.sge = 1;
%        par.run = 1;
        par.sge_queu = 'normal,bigmem';
        par.jobname  = sprintf('RS_VOI_est_double_first_%s',patient_list{subj});
        job_first_level_estimate(fspm,par);
    end
end
%%  Define contrasts

clear roi_name par
par.sge = 0;
par.run = 1;
par.sge_queu = 'normal,bigmem';
par.jobname = 'RS_VOI_con_double_first';


PositiveEffectV1= [1 -1 0 0];
PositiveEffectV2= [-1 1 0 0];

for subj = 1 : nSubj
    for ir = 1 : nROI
        fspm = gfile(outdirs{subj,ir},'SPM.mat');
        if ir == nROI
            roi_name = 'PCC';
        else
            [~,roi_name] = fileparts(fileROI{ir});
            roi_name = roi_name(1:end-10);
        end
        contrast.names={sprintf('Effect_%s_V1>V2',roi_name), sprintf('Effect_%s_V1<V2',roi_name)};
        contrast.values={PositiveEffectV1,PositiveEffectV2};
        contrast.types={'T','T'};

        par.delete_previous = 1;

        j=job_first_level_contrast(fspm,contrast,par);
    end
end
% %%script_report_edit(e);

%% Display

%%e.getOne.getModel(model_name{2}).show
