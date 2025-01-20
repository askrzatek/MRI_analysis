%% Init
clear
clc

%% contrasts definition for each model
%% Contrast : definition case 1 & 2 for both REMINARY & PARKGAME

%% Contrast : definition

LEFT_REAL_S2                            = [1 0 0 0];
LEFT_IMAGINARY_S2                       = [0 1 0 0];
RIGHT_REAL_S2                           = [0 0 1 0];
RIGHT_IMAGINARY_S2                      = [0 0 0 1];

%% T contrast

%% ANOVA2x2_LxT

contrast_T.names = {
    
    
    'LEFT_REAL_S2'
    'LEFT_IMAGINARY_S2'
    'RIGHT_REAL_S2'
    'RIGHT_IMAGINARY_S2'
    
    %'LEFT-RIGHT_S2'
    %'RIGHT-LEFT_S2'
    %'IMAGINARY-REAL_S2'
    %'REAL-IMAGINARY_S2'
    
    'POS-LATERALITY-TASK-INTERACTION_S2'
    'NEG-LATERALITY-TASK-INTERACTION_S2'
    
    %'RIGHT_REAL-LEFT_REAL_S2'
    %'RIGHT_IMAGINARY-LEFT_IMAGINARY_S2'
    
    %'LEFT_REAL-RIGHT_REAL_S2'
    %'LEFT_IMAGINARY-RIGHT_IMAGINARY_S2'
    
    'LEFT_REAL_S2 - LEFT_IMAGINARY_S2'
    'LEFT_IMAGINARY_S2 - LEFT_REAL_S2'
    
    'RIGHT_REAL_S2 - RIGHT_IMAGINARY_S2'
    'RIGHT_IMAGINARY_S2 - RIGHT_REAL_S2'
    
    
}';

contrast_T.values = {
    
    LEFT_REAL_S2        % [ 1 0 0 0]
    LEFT_IMAGINARY_S2   % [ 0 1 0 0]
    RIGHT_REAL_S2       % [ 0 0 1 0]
    RIGHT_IMAGINARY_S2  % [ 0 0 0 1] 
   
    
    %% [1 1 -1 -1]      'LEFT-RIGHT_S2'
    %(LEFT_REAL_S2 + LEFT_IMAGINARY_S2)       - (RIGHT_REAL_S2 + RIGHT_IMAGINARY_S2)
    %% [-1 -1 1 1]      'RIGHT-LEFT_S2'
    %(RIGHT_REAL_S2 + RIGHT_IMAGINARY_S2)       - (LEFT_REAL_S2 + LEFT_IMAGINARY_S2)
    %% [-1 1 -1 1]      'IMAGINARY-REAL_S2'
    %(LEFT_IMAGINARY_S2 + RIGHT_IMAGINARY_S2)  - (LEFT_REAL_S2 + RIGHT_REAL_S2)
    %% [1 -1 1 -1]      'REAL-IMAGINARY_S2'
    %(LEFT_REAL_S2 + RIGHT_REAL_S2)            - (LEFT_IMAGINARY_S2 + RIGHT_IMAGINARY_S2)
    
    % [ 1 -1 -1 1]     'POS-LATERALITY-TASK-INTERACTION_S2'
    (LEFT_REAL_S2 + RIGHT_IMAGINARY_S2)       - (RIGHT_REAL_S2 + LEFT_IMAGINARY_S2)
    % [ -1 1 1 -1]     'NEG-LATERALITY-TASK-INTERACTION_S2'
    (RIGHT_REAL_S2 + LEFT_IMAGINARY_S2)       - (LEFT_REAL_S2 + RIGHT_IMAGINARY_S2)
    
    %% [-1 0 1 0]       'RIGHT_REAL-LEFT_REAL_S2'
    %RIGHT_REAL_S2                         - LEFT_REAL_S2
    %% [ 0 -1 0 1]      'RIGHT_IMAGINARY-LEFT_IMAGINARY_S2'
    %RIGHT_IMAGINARY_S2                    - LEFT_IMAGINARY_S2
    
    %% [1 0 -1 0]       'LEFT_REAL-RIGHT_REAL_S2'
    %LEFT_REAL_S2                         - RIGHT_REAL_S2
    %% [ 0 1 0 -1]      'LEFT_IMAGINARY-RIGHT_IMAGINARY_S2'
    %LEFT_IMAGINARY_S2                    - RIGHT_IMAGINARY_S2
    
    % [-1 1 0 0]
    LEFT_IMAGINARY_S2 - LEFT_REAL_S2

    % [1 -1 0 0]
    LEFT_REAL_S2 - LEFT_IMAGINARY_S2
    
    % [0 0 1 -1]
    RIGHT_REAL_S2 - RIGHT_IMAGINARY_S2
    
    % [0 0 -1 1]
    RIGHT_IMAGINARY_S2 - RIGHT_REAL_S2
    
}';


contrast_T.types = cat(1,repmat({'T'},[1 length(contrast_T.names)]));


%% F contrast

contrast_F.names = {
    
    'Avg_Condition_Effect_S2'
    'Main_LATERALITY_Effect_LEFT_S2'
    'Main_LATERALITY_Effect_RIGHT_S2'
    'Main_TASK_Effect_IMA_S2'
    'Main_TASK_Effect_REAL_S2'
    'Main_LxT_INTERACTION_Effect_S2'
    
    'Main_LEFT_REAL_S2'
    'Main_LEFT_IMAGINARY_S2'
    'Main_RIGHT_REAL_S2'
    'Main_RIGHT_IMAGINARY_S2'
    
}';
%%% CHECK THAT !!! EQUILIBRIUM aspect !!! should it not be 2* (LEFT_REAL_S2 + LEFT_IMAGINARY_S2) ???

contrast_F.values = {
    %Avg_Condition_Effect_S2            = [1 1 1 1];
    LEFT_REAL_S2 + LEFT_IMAGINARY_S2 + RIGHT_REAL_S2 + RIGHT_IMAGINARY_S2
    %Main_LATERALITY_Effect_LEFT_S2             = [1 1 0 0];
    LEFT_REAL_S2                 + LEFT_IMAGINARY_S2
    %Main_LATERALITY_Effect_RIGHT_S2             = [0 0 1 1];
    RIGHT_REAL_S2                 + RIGHT_IMAGINARY_S2
    %Main_TASK_Effect_IMA_S2                = [0 1 0 1];
    LEFT_IMAGINARY_S2            + RIGHT_IMAGINARY_S2
    %Main_TASK_Effect_REAL_S2                = [1 0 1 0];
    LEFT_REAL_S2            + RIGHT_REAL_S2
    %Main_SxT_INTERACTION_Effect_S2     = [1 0 0 1];
    LEFT_REAL_S2                 + RIGHT_IMAGINARY_S2
   
    %'Main_LEFT_REAL_S2'            = [1 0 0 0]
    LEFT_REAL_S2
    %'Main_LEFT_IMAGINARY_S2'       = [0 1 0 0]
    LEFT_IMAGINARY_S2
    %'Main_RIGHT_REAL_S2'           = [0 0 1 0]
    RIGHT_REAL_S2
    %'Main_RIGHT_IMAGINARY_S2'      = [0 0 0 1]
    RIGHT_IMAGINARY_S2
    
    
}';


contrast_F.types = cat(1,repmat({'F'},[1 length(contrast_F.names)]));


contrast_2x2_LxT.names  = [contrast_F.names  contrast_T.names];
contrast_2x2_LxT.values = [contrast_F.values contrast_T.values];
contrast_2x2_LxT.types  = [contrast_F.types  contrast_T.types];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save the model_names to create directories and their corresponding contrast structures

model_name     = {'ANOVA2x2_LxT'};
model_contrast = {contrast_2x2_LxT};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fetch dirs for SPM
%% Load files from multiple folders
main_dir = fullfile('/network/iss/cenir/analyse/irm/users/anna.skrzatek','nifti');
cd (main_dir)
% e_PARKGAME = exam(main_dir,'PARKGAME');
% e_REMINARY = exam(main_dir,'REMINARY_\w{2}_');

%e_REMINARY_S2 = exam(main_dir,'REMINARY_\w{2}_.*2$');
%e = {e_REMINARY_S2};
%dirstat = r_mkdir(main_dir, 'secondlevel_ACTIVATION_REM_S2');
%dirgroup = dirstat;

e_PARKGAME_S2_a = exam(main_dir,'PARKGAME.*2_a'); % taking into account all S2 patients
e_PARKGAME_S2_c = exam(main_dir,'PARKGAME.*2_c'); % taking into account all S2 patients
e = {e_PARKGAME_S2_a, e_PARKGAME_S2_c};
dirstat = r_mkdir(main_dir, 'secondlevel_ACTIVATION_PARK_S2');

dirgroup = {'PARK_a', 'PARK_c'};
%dirgroup = r_mkdir(char(dirstat), {'PARKGAME_S2'});

done = 0;
%done = job_con_smooth('s',4); % comment this line if you don't want to smooth your contrast data
if done ==1
    for i = 1:length(e)
        %e{i}.explore
        %'REMINARY_\w{2}_.*1$'
        e{i}.addSerie('model_tedana$','contrasts',1)

        e{i}.getSerie('contrasts').addVolume('^scon_0008','REAL_L',1)
        e{i}.getSerie('contrasts').addVolume('^scon_0010','IMA_L',1)
        e{i}.getSerie('contrasts').addVolume('^scon_0009','REAL_R',1)
        e{i}.getSerie('contrasts').addVolume('^scon_0011','IMA_R',1)

        [ec_1st, ei_1st] = e{i}.removeIncomplete;
        e{i} = ec_1st;

        e{i}.explore
        %e{2}.explore
       dirfig = 'auto_figures_smooth';
    end
else
    for i = 1:length(e)
        %e{i}.explore
        %'REMINARY_\w{2}_.*1$'
        e{i}.addSerie('model_tedana$','contrasts',1)

        e{i}.getSerie('contrasts').addVolume('^con_0008','REAL_L',1)
        e{i}.getSerie('contrasts').addVolume('^con_0010','IMA_L',1)
        e{i}.getSerie('contrasts').addVolume('^con_0009','REAL_R',1)
        e{i}.getSerie('contrasts').addVolume('^con_0011','IMA_R',1)

        [ec_1st, ei_1st] = e{i}.removeIncomplete;
        e{i} = ec_1st;

        e{i}.explore
        %e{2}.explore
        dirfig = 'auto_figures';
    end
end


par.fake = 0;
par.redo = 0;
par.verbose = 2;

%% create a result directories
%dirout = r_mkdir(dirstat{1}, model_name);

%% Probably a MODEL LOOP will start here with model_name and model_contrast varying with iterations
for group=1:length(dirgroup)
    for imod=1:length(model_name)
        if group == 1 % keep if we still need 2 exams per group
            dirout = r_mkdir(dirgroup{group}, model_name);
        else
            dirout = r_mkdir(dirgroup{2}, model_name); % works only if we have no more than 2 groups ergo 2 protocols
        end
        model_dir = cellstr(dirout{imod});
        fact1 = 'Hand';
        fact2 = 'Task';
        l11 = e{group}.getSerie('contrasts').getVolume('REAL_L').toJob;
        l12 = e{group}.getSerie('contrasts').getVolume('IMA_L').toJob; 
        l21 = e{group}.getSerie('contrasts').getVolume('REAL_R').toJob;
        l22 = e{group}.getSerie('contrasts').getVolume('IMA_R').toJob;
        addpath '/network/iss/cenir/analyse/irm/users/anna.skrzatek/'
%% Job define model

        par.fake = 0;
        par.redo = 0;
        par.verbose = 2;
        par.run = 1;
%   par.file_reg = '^con'; % it doesn't mean anything - you can put whatever you want, it's l11, l12 etc that define the files
        par.mask_thr = 0.07;

        job_second_level_specify(fact1,fact2,l11,l12,l21,l22,model_dir,par);

%% Estimate
        fspm = fullfile(model_dir,'SPM.mat');
        job_second_level_estimate(fspm);

%% Contrast : write

        par.run = 1;
        par.display = 0;

    % par.sessrep = 'both';
        par.sessrep = 'none';

        par.delete_previous = 1;
    %par.report          = 1; % error on non designating spm mode
    %(?)
        par.report          = 0;
            
        job_second_level_contrast(fspm,model_contrast{imod},par);
%% create a folder for figures before creating figures with MRIcroGL
        mkdir(dirout{imod},dirfig);
    end
end
%save e

%% Display

%!linux command for mricrogl script
%!/network/iss/cenir/software/irm/mricrogl_lx/MRIcroGL '/home/anna.skrzatek/data/nifti/second_level_p001_auto_less.gls'

%addpath /network/iss/cenir/software/irm/spm12/toolbox/marsbar/


