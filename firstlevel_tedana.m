%% Init

clear
clc

addpath /home/anna.skrzatek/data/StimTemplate
cd      /home/anna.skrzatek/data

%main_dir = fullfile(pwd,'nifti');

% stim_dir = fullfile(pwd,'/nifti_test/behav_test');
% main_dir = fullfile(pwd,'/nifti_test/ben');

stim_dir = fullfile(pwd,'/behav');
main_dir = fullfile(pwd,'/nifti_test');

%model_name = 'model_meica';
%model_name = 'model_tedana';
%model_name = {'model_tedana', 'model_ts_tapas'};
%model_name = {'smodel_tedana','smodel_ts_tapas'};%, 'smodel_dn_tapas'};
model_name  = {'rsmodel_ts_tapas'};

%% fetch dirs
cd (main_dir)
load e % why doesn't the e-object work?!

%e = e(2:length(e));

dir_func  = e.getSerie('run_ACTIVATION') .toJob;
dir_anat = e.getSerie('anat').toJob(0); % useful only if individual display in the end (I believe)

par.fake = 0;
par.redo = 0;
par.verbose = 2;

%% add the if condition here when we only need tedana results

%model_outdir1 = e.mkdir(model_name);
model_outdir1 = e.mkdir(model_name{1});
%model_outdir2 = e.mkdir(model_name{2});
%model_outdir3 = e.mkdir(model_name{3});


%% Fetch onset & mask

e.getSerie('run_ACTIVATION').addStim(stim_dir, 'MRI_run\d{2}_SPM.mat', 'run', 1 )
%[ec_stim, ei_stim] = e.removeIncomplete;
stim_files = e.getSerie('run_ACTIVATION').getStim.toJob(0);
e.getSerie('run_ACTIVATION').addVolume('^w.*mask.nii','mask',1)

% e.explore

%% Make symbolic links from tedana_vtd_mle dir to run dir based on job_meica_afni symbolic link creation 
%par.subdir        = 'tedana_vtd_mle';
par.subdir        = 'tedana009a1_vtd';
par.sge = 0;
par.run = 1;

% par.warp_file_reg = '^wdn';
% job_symbolic_child_to_parent(dir_func, par);
% 
% par.warp_file_reg = '^s5wts';
% job_symbolic_child_to_parent(dir_func, par);
% 
% par.warp_file_reg = '^s5wdn';
% job_symbolic_child_to_parent(dir_func, par);

par.warp_file_reg = '^s6wts';
job_symbolic_child_to_parent(dir_func, par);

% par.warp_file_reg = '^s6wdn';
% job_symbolic_child_to_parent(dir_func, par);

%% make a symbolic link of rp_spm.txt and of multiple_regressors to dir_func


for i= 1:length(e)
    par.subdir = 'wts';
    par.regfile_regex = 'multiple_regressors.txt';
    wd = e(i).getSerie('run_ACTIVATION').path;
    reg_dir = char(get_subdir_regex(wd, par.subdir));
    A_src = fullfile(reg_dir, par.regfile_regex);
    regfile_out = sprintf('%s_%s',par.subdir,par.regfile_regex);
    A_dst = fullfile(wd, regfile_out);
    
    par.redo = 0;
    par.verbose = 2;
    par.run = 0;
    par.run = 1;
    par.jobname = sprintf('job_symbolic_link');
    %par.jobname = sprintf('%s_%s_%s', 'job_symbolic_link', wd(end-46:end-16), par.subdir);

    %% reste à tester    
    if ~exist(fullfile(wd, regfile_out))
        [job_session(i)] = r_movefile(A_src, A_dst, 'linkn', par);
        job = [job_session];
    else
        clear par
    end
end
%job = do_cmd_sge(job, par);


%% Job define model
clear par
par.sge = 0;
par.run = 0;
par.display = 0;
%par.mask = e.getSerie('run_ACTIVATION').getVolume('mask').toJob;

par.mask_thr = 0.1;
%par.mask_thr = 0.07;
par.TR = 1.6;
par.rp = 1;

%% Classic tedana

addpath /home/anna.skrzatek/matvol/SPM/firstlevel/

% par.file_reg = '^s6wdn.*nii';
% par.rp_regex = 'rp.*txt';
% job1 = job_first_level_specify_xnat(dir_func, model_outdir1, stim_files, par);

%% ts TAPAS
%par.rp = 1;
par.file_reg = '^s6wts.*nii';
par.rp_regex = 'wts_multiple_regressors.txt';
job1 = job_first_level_specify_xnat(dir_func, model_outdir1, stim_files, par);

% %% tedana + TAPAS
% 
% par.file_reg = '^s6wdn.*nii';
% 
% par.rp_regex = 'wdn_multiple_regressors.txt';
% job3 = job_first_level_specify(dir_func, model_outdir3, stim_files, par);

jobs = job1 ; %job2 ]; %job3 ];

par.sge = 1;
par.run = 0;
par.display = 0;
par.jobname = 'spm_single_glm_def';

job_ending_rountines(jobs,[],par);

return

%% Estimate

e.addModel(model_name{1},model_name{1});
%e.addModel(model_name{2},model_name{2});
%e.addModel(model_name{3},model_name{3});
%e.addModel('^model_tedana','model_tedana');

%ei = e(34);
%e = e(1:33);
%save ('e','e')
clear par

par.run = 0;
par.sge = 1;
par.display = 0;
par.jobname = 'spm_single_glm_est'
%fspm = e.getModel(cellstr2regex({model_name{1}, model_name{2}},1)).removeEmpty.toJob;
fspm = e.getModel(model_name{1}).removeEmpty.toJob; %% à tester

job_first_level_estimate(fspm,par);


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

%% complementary individual contrasts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'REAL_Left - REAL_Right'
% 'REAL_Right - REAL_Left'
% 'IMAGINARY_Left - IMAGINARY_Right'
% 'IMAGINARY_Right - IMAGINARY_Left'
% 
% 'REAL_Left - IMAGINARY_Left'
% 'IMAGINARY_Left - REAL_Left'
% 'REAL_Right - IMAGINARY_Right'
% 'IMAGINARY_Right - REAL_Right'
% 
% 'IMAGINARY total - REAL total'
% 'REAL total       - IMAGINARY total'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 'Instruction - Rest'
% 'Instruction - Total Activation'
% 'Total Activation - Instruction'
% 'Total_REAL_Activation - Instruction'
% 'REAL_Left - Instruction'
% 'REAL_Right - Instruction'
% 'Total_IMAGINARY_Activation - Instruction'
% 'IMAGINARY_Left - Instruction'
% 'IMAGINARY_Right - Instruction'

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

%% complementary individual contrasts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REAL_Left       - REAL_Right
% REAL_Right      - REAL_Left
% IMAGINARY_Left  - IMAGINARY_Right
% IMAGINARY_Right - IMAGINARY_Left
% 
% REAL_Left       - IMAGINARY_Left
% IMAGINARY_Left  - REAL_Left
% REAL_Right      - IMAGINARY_Right
% IMAGINARY_Right - REAL_Right
% (IMAGINARY_Left + IMAGINARY_Right) - (REAL_Left + REAL_Right)
% (REAL_Left + REAL_Right)           - (IMAGINARY_Left + IMAGINARY_Right)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Instruction     - Rest
% 4*Instruction   - (REAL_Left + REAL_Right + IMAGINARY_Left + IMAGINARY_Right)
% (REAL_Left + REAL_Right + IMAGINARY_Left + IMAGINARY_Right) - 4*Instruction
% (REAL_Left + REAL_Right)                                    - 2*Instruction
% REAL_Left       - Instruction
% REAL_Right      - Instruction
% (IMAGINARY_Left + IMAGINARY_Right)                          - 2*Instruction
% IMAGINARY_Left  - Instruction
% IMAGINARY_Right - Instruction

}';


contrast_T.types = cat(1,repmat({'T'},[1 length(contrast_T.names)]));


% F contrast

contrast_F.names = {'F-all'}';

contrast_F.values = {eye(6)}';


contrast_F.types = cat(1,repmat({'F'},[1 length(contrast_F.names)]));


contrast.names  = [contrast_F.names  contrast_T.names];
contrast.values = [contrast_F.values contrast_T.values];
contrast.types  = [contrast_F.types  contrast_T.types];

% par.contrast = contrast;
% par.jobname  = 'spm_new_syntax_it_would_be_to_easy_otherwise';
% job_contrast(model_outdir1, par)

%% Contrast : write
clear par

par.sge = 0;
par.run = 1;
par.display = 0;
par.jobname = 'spm_single_glm_con';

% par.sessrep = 'both';
par.sessrep = 'none';

par.delete_previous = 1;
par.report          = 0;

%fspm = e.getModel(cellstr2regex({model_name{1}, model_name{2}},1)).removeEmpty.toJob;
fspm = e.getModel(model_name{1}).removeEmpty.toJob; % à tester

job_first_level_contrast(fspm,contrast,par);

%% add the model and the contrasts to the e object

%e.addSerie('^model','model',3);
e.getSerie('model').addVolume('^SPM','m',1);
e.getSerie('model').addVolume('^spm[TF]_\d{4}','con',length(contrast.names));

% %% Display and save results for all subjects, for one contrast example
% 
% e.getSerie('anat').addVolume('^wp0','t1',1);
% t1 = e.getSerie('anat').getVolume('t1');
% Coordlist.values = cat(2,[0.0; 0.0; 0.0], [34.0; -24.0; 68.0], [-34.0; -24.0; 68.0], [-4.0; -48.0; -24.0], [0.0; -24.0; 10.0]); 
% Coordlist.names = {'centre'; 'RIGHT SM'; 'LEFT SM'; 'CEREBELLUM'; 'MOTOR BI'};
% %mkdir('/home/anna.skrzatek/Desktop','auto_figures_first_level_con4');
% output_dir = '/home/anna.skrzatek/Desktop/auto_figures_first_level_con4/';
% wd = '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/';
% 
% for n =1:length(fspm)
%     mdir = e(n).getSerie('model').path;
%     cd (mdir)
%     load SPM.mat
%     job_spm_single_results_display({fspm{n}}, 4, t1(n).path, Coordlist,  output_dir, SPM, wd)
% end
% % 

%% Save the new e object

%cd /network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/
%save e

%% Create figures
clear par

fsub = e.gpath;
Coordlist.values = cat(2,[0.0; 0.0; 0.0], [34.0; -24.0; 68.0], [-34.0; -24.0; 68.0], [-4.0; -48.0; -24.0], [0.0; -24.0; 10.0]); % EXAMPLE
Coordlist.names = {'centre'; 'RIGHT SM'; 'LEFT SM'; 'CEREBELLUM'; 'MOTOR BI'}; % EXAMPLE

par.subdir       = 'smodel_tedana';
par.anat_dir_reg = 'S\d{2}_t1mpr_S256_0_8iso_p2$';
par.output_dir   = 'smooth6_tedana_figures';
par.extent       = 5;
par.thresh       = 0.001;
par.fixscale     = 1;
par.minscale     = 0;
par.maxscale     = 20;

%conlist = {'F-all','REAL_Left','REAL_Right','IMAGINARY_Left','IMAGINARY_Right','REAL_Left - Rest','REAL_Right - Rest','IMAGINARY_Left - Rest','IMAGINARY_Right - Rest'};
conlist = {'REAL_Right - Rest','IMAGINARY_Right - Rest'};
addpath /home/anna.skrzatek/matvol/SPM/

for icon = 1 : length(conlist)
   par.conname  = conlist{icon};
   job_spm_single_results_display(fsub, Coordlist, par);
end

close all

% %%script_report_edit(e);

%% Display

%%e.getOne.getModel(model_name{2}).show
