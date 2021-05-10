%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% script for creating the regression models between fMRI & other statistical measures
%% ex. Clinical, Gait or Game Data regression for each fMRI condition
%% Anna SKRZATEK May 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all

%% Initialise

main_dir = fullfile('/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek','/nifti_test');
cd (main_dir)

ACTION = 1;
RS = 0;

%% define or create output directory
% %RegDir = mkdir('resliced_multiple_regression')
% RegDir = fullfile(main_dir, '/resliced_ACT_clinic');
RegDir = fullfile(main_dir, '/resliced_ACT_gait');

%% Regressors definition
Clinic = {'AXIAL','GABS','UPDRSIII','UPDRSIII-AXIAL'};
Gait  = {'Rapid/rAPA_AP_duration','Rapid/rDA_duration','Rapid/rStepSize','Spontaneous/sAPA_AP_duration','Spontaneous/sDA_duration','Spontaneous/sStepSize'};
Gait_up  = {'Rapid','Spontaneous'};
%Gait_sub  = {'rAPA_AP_duration','Rapid/rDA_duration','Rapid/rStepSize','Spontaneous/sAPA_AP_duration','Spontaneous/sDA_duration','Spontaneous/sStepSize'};


model = Gait_up;
for ivar = 1 : length(model)
%     mkdir(RegDir,Clinic{ivar});
%     StatDir{ivar} = fullfile(RegDir,Clinic{ivar});
    mkdir(RegDir,Gait{ivar});
    StatDir{ivar} = fullfile(RegDir,Gait{ivar});
    StatObj{ivar} = exam(RegDir,Gait_up{ivar},'^[r,s]');
end

Stat = StatObj{1} + StatObj{2};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ACTION
%% Choice of fMRI conditions //scans, condition names 
if ACTION
    Conditions = {'IL','IR','RL','RR'};
    Sessions = {'V1','V2','V2-V1'};
    for iC = 1 : length(Conditions)
        for iS = 1 : length(Sessions)
            mkdir(StatDir{:},char(sprintf('%s_%s',Conditions{iC},Sessions{iS})))
            Stat.mkdir(sprintf('%s_%s',Conditions{iC},Sessions{iS})) .toJob;
            
        end
    end

%% Models specification

%%
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RS
%% Choice of fMRI conditions //scans, condition names 
if RS
    

%% Models specification

%%
 
end

%% Models estimation

fspm = addsuffixtofilenames(StatDir, 'SPM.mat');

clear par
%par.run = 1;
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.jobname  = 'spm_first_level_est_RS_wbet';
job_first_level_estimate(fspm,par)