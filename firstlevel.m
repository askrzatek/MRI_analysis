%% Init

clear
clc

addpath /home/anna.skrzatek/data/StimTemplate

%load e % if in the correct directory - cd wd before

main_dir = fullfile(pwd,'nifti');
stim_dir = fullfile(pwd,'behav');

model_name = 'model_meica';
%model_name = 'model_tedana';

%% fetch dirs

dir_func  = e.getSerie('run_ACTIVATION') .toJob;
dir_anat = e.getSerie('anat').toJob(0);

par.fake = 0;
par.redo = 0;
par.verbose = 2;

model_dir = e.mkdir(model_name);


%% Fetch onset

e.getSerie('run_ACTIVATION').addStim(stim_dir, 'MRI_run\d{2}_SPM.mat', 'run', 1 )
stim_files = e.getSerie('run_ACTIVATION').getStim.toJob(0);

% e.explore


%% Job define model

par.run = 1;
%par.pct = 1;
% par.TR = 1.6;

par.file_reg = '^wdn';

% par.rp = 1;
par.mask_thr = 0.07;

%% Define model

job_first_level_specify(dir_func,model_dir,stim_files,par)


%% Estimate

fspm = e.addModel(model_name,model_name);
job_first_level_estimate(fspm,par)


%% Contrast : definition

Rest            = [1 0 0 0 0 0];
REAL_Left       = [0 1 0 0 0 0];
REAL_Right      = [0 0 1 0 0 0];
IMAGINARY_Left  = [0 0 0 1 0 0];
IMAGINARY_Right = [0 0 0 0 1 0];
Instruction     = [0 0 0 0 0 1];

% T contrast

contrast_T.names = {

'Rest'
'REAL_Left'
'REAL_Right'
'IMAGINARY_Left'
'IMAGINARY_Right'
'Instruction'

'REAL_Left - Rest'
'REAL_Right - Rest'
'IMAGINARY_Left - Rest'
'IMAGINARY_Right - Rest'

%'REAL_Left - REAL_Right'
%'REAL_Right - REAL_Left'
%'IMAGINARY_Left - IMAGINARY_Right'
%'IMAGINARY_Right - IMAGINARY_Left'

%'REAL_Left - IMAGINARY_Left'
%'IMAGINARY_Left - REAL_Left'
%'REAL_Right - IMAGINARY_Right'
%'IMAGINARY_Right - REAL_Right'

'IMAGINARY total - REAL total'
'REAL total       - IMAGINARY total'


%'Instruction - Rest'
%'Instruction - Total Activation'
%'Total Activation - Instruction'
%'Total_REAL_Activation - Instruction'
%'REAL_Left - Instruction'
%'REAL_Right - Instruction'
%'Total_IMAGINARY_Activation - Instruction'
%'IMAGINARY_Left - Instruction'
%'IMAGINARY_Right - Instruction'

}';

contrast_T.values = {
    
Rest
REAL_Left
REAL_Right
IMAGINARY_Left
IMAGINARY_Right
Instruction

REAL_Left       - Rest
REAL_Right      - Rest
IMAGINARY_Left  - Rest
IMAGINARY_Right - Rest

%REAL_Left       - REAL_Right
%REAL_Right      - REAL_Left
%IMAGINARY_Left  - IMAGINARY_Right
%IMAGINARY_Right - IMAGINARY_Left

%REAL_Left       - IMAGINARY_Left
%IMAGINARY_Left  - REAL_Left
%REAL_Right      - IMAGINARY_Right
%IMAGINARY_Right - REAL_Right

(IMAGINARY_Left + IMAGINARY_Right) - (REAL_Left + REAL_Right)
(REAL_Left + REAL_Right)           - (IMAGINARY_Left + IMAGINARY_Right)

%Instruction     - Rest
%4*Instruction   - (REAL_Left + REAL_Right + IMAGINARY_Left + IMAGINARY_Right)
%(REAL_Left + REAL_Right + IMAGINARY_Left + IMAGINARY_Right) - 4*Instruction
%(REAL_Left + REAL_Right)                                    - 2*Instruction
%REAL_Left       - Instruction
%REAL_Right      - Instruction
%(IMAGINARY_Left + IMAGINARY_Right)                          - 2*Instruction
%IMAGINARY_Left  - Instruction
%IMAGINARY_Right - Instruction

}';


contrast_T.types = cat(1,repmat({'T'},[1 length(contrast_T.names)]));


% F contrast

contrast_F.names = {

'F-all'

}';

contrast_F.values = {
    
eye(6)

}';


contrast_F.types = cat(1,repmat({'F'},[1 length(contrast_F.names)]));


contrast.names  = [contrast_F.names  contrast_T.names];
contrast.values = [contrast_F.values contrast_T.values];
contrast.types  = [contrast_F.types  contrast_T.types];


%% Contrast : write

par.run = 1;
par.display = 0;

% par.sessrep = 'both';
par.sessrep = 'none';

par.delete_previous = 1;
par.report          = 0;
job_first_level_contrast(fspm,contrast,par);

% add the model and the contrasts to the e object

e.addSerie('^model','model',1);
e.getSerie('model').addVolume('^SPM','m',1);
e.getSerie('model').addVolume('^spm[TF]_\d{4}','con',length(contrast.names));

% %% Display and save results for all subjects, for one contrast example
% 
% e.getSerie('anat').addVolume('^wms','t1',1);
% t1 = e.getSerie('anat').getVolume('t1');
% Coordlist.values = cat(2,[0.0; 0.0; 0.0], [34.0; -24.0; 68.0], [-34.0; -24.0; 68.0], [-4.0; -48.0; -24.0], [0.0; -24.0; 10.0]); 
% Coordlist.names = {'centre'; 'RIGHT SM'; 'LEFT SM'; 'CEREBELLUM'; 'MOTOR BI'};
% mkdir('/home/anna.skrzatek/Desktop','auto_figures_first_level_con4');
% output_dir = '/home/anna.skrzatek/Desktop/auto_figures_first_level_con4/';
% wd = '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/';
% 
% for n =1:length(fspm)
%     mdir = e(n).getSerie('model').path;
%     cd (mdir)
%     load SPM.mat
%     job_spm_single_results_display({fspm{n}}, 4, t1(n).path, Coordlist,  output_dir, SPM, wd)
% end
% % 
cd /network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/
save e


%script_report_edit(e);

%% Display

%e.getOne.getModel(model_name).show
