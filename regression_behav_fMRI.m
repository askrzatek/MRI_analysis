%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% script for creating the regression models between fMRI & other statistical measures
%% ex. Clinical, Gait or Game Data regression for each fMRI condition
%% Anna SKRZATEK May 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
par.jobname = 'regression_model_spec'

%% Initialise

main_dir = fullfile('/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek','/nifti_test');
cd (main_dir)


ACTION = 1;
RS = 0;

% define input directory
InputfMRI = fullfile(main_dir,'/double_run_resliced');
InputRS   = fullfile(main_dir,'RS');

%% define or create output directory
% %RegDir = mkdir('resliced_multiple_regression')

% RegDir = fullfile(main_dir, '/resliced_ACT_clinic');
RegDir = fullfile(main_dir, '/resliced_ACT_gait');

%% Regressors definition
%% AGE
covars{1} = [70
             74
             64
             76
             61
             75
             66
             72
             68
             68
             72
             56
             57];
%% GENDER
covars{2} =[1
            1
            1
            2
            1
            2
            2
            1
            2
            2
            2
            2
            1];

%% Target regressors definition

Clinic = {'AXIAL','GABS','UPDRSIII','UPDRSIII-AXIAL'};
Gait  = {'Rapid/rAPA_AP_duration','Rapid/rDA_duration','Rapid/rStepSize','Spontaneous/sAPA_AP_duration','Spontaneous/sDA_duration','Spontaneous/sStepSize'};
Gait_up  = {'Rapid','Spontaneous'};


model = Gait_up;
for ivar = 1 : length(model)
    mkdir(RegDir,sprintf('%s',model{ivar}); % universal creation of directories depending on the chosen model
    % StatDir{ivar} = fullfile(RegDir,Gait{ivar});
    StatObj{ivar} = exam(RegDir,sprintf('%s$',Clinic{ivar}); % we have two of them containing UPDRSIII
    StatObj{ivar} = exam(RegDir,Gait_up{ivar},'^r'); % the Spontaneous gait condition is one subject shorter - needs a separate processing
end

Stat = sum(StatObj)
Stat.explore


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ACTION
%% Choice of fMRI conditions //scans, condition names 
if ACTION
    Conditions  = {'IL','IR','RL','RR'};
    Sessions    = {'V1','V2','V2-V1'};
% input different groups distinction and paths for the cons we need
    fMRIObj = exam(InputfMRI,'firstlevel') % check directory names
    fMRIObj.addSerie('PARK.*a$','a',8)
    fMRIObj.addSerie('PARK.*c$','c',6)
% add volumes for each contrast - each individual

end

if RS
    ROIs        = {'PCC','Motor_L','Motor_R','SMA_L','SMA_R'}; 
    Conditions  = ROIs;
    Sessions    = {'V1','V2','V2-V1'};
% input different groups distinction and paths for the cons we need
    RSObj = exam(InputRS,'firstlevel') % check directory names
    RSObj.addSerie('PARK.*a$','a',8)
    RSObj.addSerie('PARK.*c$','c',6)
% add volumes for each contrast - each individual

end

    %% this part will probably become in common with the RS processing

    for iC = 1 : length(Conditions)
        for iS = 1 : length(Sessions)
            %mkdir(StatDir{:},char(sprintf('%s_%s',Conditions{iC},Sessions{iS})))
            Stat.mkdir(sprintf('%s_%s',Conditions{iC},Sessions{iS}))
            Stat.addSerie(sprintf('%s_%s',Conditions{iC},Sessions{iS}),sprintf('%s_%s',Conditions{iC},Sessions{iS}))
        end
    end

    outdirs = Stat.getSerie(sprintf('%s_%s',Conditions{iC},Sessions{iS})) .toJob;


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