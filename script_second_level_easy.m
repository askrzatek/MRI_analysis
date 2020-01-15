clear all
close all
%%
main_dir = fullfile('/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek','nifti');

e_REMINARYS1 = exam(main_dir,'REMINARY_\w{2}_.*1$');
e_REMINARYS2 = exam(main_dir,'REMINARY_\w{2}_.*2$');
e_PARKGAMES2 = exam(main_dir,'PARKGAME.*2$');
e_PARKGAMES1 = exam(main_dir,'PARKGAME.*1$');

%% case 1 REMINARY_ANOVA_2x2
fact1 = 'SESSION';
fact2 = 'TASK';

e_REMINARYS1.addSerie('model_meica','source',1);
l11 = e_REMINARYS1.getSerie('source').addVolume('^con_\d{2}08','realRightS1',1);

e_REMINARYS2.addSerie('model_meica','source',1);
l21 = e_REMINARYS2.getSerie('source').addVolume('^con_\d{2}08','realRightS2',1);

e_REMINARYS1.addSerie('model_meica','source',1);
l12 = e_REMINARYS1.getSerie('source').addVolume('^con_\d{2}10','imaRightS1',1);

e_REMINARYS2.addSerie('model_meica','source',1);
l22 = e_REMINARYS2.getSerie('source').addVolume('^con_\d{2}10','imaRightS2',1);

mkdir(main_dir,'secondlevel_ACTIVATION');
dirout = get_subdir_regex(main_dir,'^second.*N$');

%% second level specify
par.fake = 0;
par.redo = 0;
par.verbose = 2;
par.run = 1;
par.file_reg = '^con';
par.mask_thr = 0.07;

job_second_level_specify(fact1,fact2,l11,l12,l21,l22,dirout,par);

%% model estimate
fspm = fullfile(dirout,'SPM.mat');
job_second_level_estimate(fspm);

%% Contrast : definition

S1REAL_R                           = [1 0 0 0];
S1IMAGINARY_R                      = [0 1 0 0];
S2REAL_R                           = [0 0 1 0];
S2IMAGINARY_R                      = [0 0 0 1];

S1REAL_L                           = [1 0 0 0];
S1IMAGINARY_L                      = [0 1 0 0];
S2REAL_L                           = [0 0 1 0];
S2IMAGINARY_L                      = [0 0 0 1];

%% T contrast
% ANOVA2x2_SESSION

contrast_T.names = {
    
    'SESSION1-SESSION2_R'
    'SESSION2-SESSION1_R'
    'IMAGINARY-REAL_R'
    'REAL-IMAGINARY_R'
    
    'POS-SESSION-TASK-INTERACTION_R'
    'NEG-SESSION-TASK-INTERACTION_R'
    
    'S2REAL-S1REAL_R'
    'S2IMAGINARY-S1IMAGINARY_R'
    
    'S1REAL-S2REAL_R'
    'S1IMAGINARY-S2IMAGINARY_R'

}';

contrast_T.values = {
    
% [1 1 -1 -1]
(S1REAL_R + S1IMAGINARY_R)       - (S2REAL_R + S2IMAGINARY_R)
% [-1 -1 1 1]
(S2REAL_R + S2IMAGINARY_R)       - (S1REAL_R + S1IMAGINARY_R)
% [-1 1 -1 1]
(S1IMAGINARY_R + S2IMAGINARY_R)  - (S1REAL_R + S2REAL_R)
% [1 -1 1 -1]
(S1REAL_R + S2REAL_R)            - (S1IMAGINARY_R + S2IMAGINARY_R)
    
% [1 -1 -1 1]
(S1REAL_R + S2IMAGINARY_R)       - (S2REAL_R + S1IMAGINARY_R)
% [ -1 1 1 -1]
(S2REAL_R + S1IMAGINARY_R)       - (S1REAL_R + S2IMAGINARY_R)

% [-1 0 1 0]
S2REAL_R                         - S1REAL_R
% [ 0 -1 0 1]
S2IMAGINARY_R                    - S1IMAGINARY_R

% [1 0 -1 0]
S1REAL_R                         - S2REAL_R
% [0 1 0 -1]
S1IMAGINARY_R                    - S2IMAGINARY_R

}';


contrast_T.types = cat(1,repmat({'T'},[1 length(contrast_T.names)]));


%% F contrast

contrast_F.names = {
    
    'Avg_Condition_Effect_R'
    'Main_SESSION_Effect_R'
    'Main_TASK_Effect_R'
    'Main_SxT_INTERACTION_Effect_R'
    
}';
%%% CHECK THAT !!! EQUILIBRIUM aspect !!! should it not be 2* (S1REAL_R + S1IMAGINARY_R) ???

contrast_F.values = {
%Avg_Condition_Effect_R            = [1 1 1 1];
    S1REAL_R + S1IMAGINARY_R + S2REAL_R + S2IMAGINARY_R
%Main_SESSION_Effect_R             = [1 1 0 0];
    S1REAL_R                 + S1IMAGINARY_R
%Main_TASK_Effect_R                = [0 1 0 1];
    S1IMAGINARY_R            + S2IMAGINARY_R
%Main_SxT_INTERACTION_Effect_R     = [1 0 0 1];
    S1REAL_R                 + S2IMAGINARY_R

}';


contrast_F.types = cat(1,repmat({'F'},[1 length(contrast_F.names)]));


contrast_2x2Session.names  = [contrast_F.names  contrast_T.names];
contrast_2x2Session.values = [contrast_F.values contrast_T.values];
contrast_2x2Session.types  = [contrast_F.types  contrast_T.types];

%%
par.run = 1;
par.display = 0;

% par.sessrep = 'both';
par.sessrep = 'none';

par.delete_previous = 1;
par.report          = 0;

job_second_level_contrast(fspm,contrast_2x2Session,par);
