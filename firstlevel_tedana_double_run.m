%% Init

clear
clc

addpath /network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/StimTemplate
addpath /home/anna.skrzatek/matvol/SPM/firstlevel/
addpath('/home/anna.skrzatek/MRI_analysis/')


cd      /network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/

%main_dir = fullfile(pwd,'nifti');
% 
% stim_dir = fullfile(pwd,'/nifti_test/behav_test');
% main_dir = fullfile(pwd,'/nifti_test/ben');

%stim_dir = fullfile(pwd,'/behav');
main_dir = fullfile(pwd,'/nifti_test');
stim_dir = fullfile(main_dir,'electrophy');

%model_name = 'model_meica';
%model_name = 'model_tedana';
%model_name = {'model_tedana', 'model_ts_tapas'};
model_name = {'full_sts_tapas_doublerun_resliced'};%, 'smodel_ts_tapas', 'smodel_dn_tapas'};

%% fetch dirs
cd (main_dir)

patient_regex = {'PARKGAMEII.*NB.*_a$','PARKGAMEII.*BM.*_a$','PARKGAMEII.*SM.*_c$','PARKGAMEII.*SD.*_a$','PARKGAMEII.*JR.*_a$','PARKGAMEII.*LJ.*_c$','PARKGAMEII.*CA.*_a$','PARKGAMEII.*PC.*_c$','PARKGAMEII.*DD.*','PARKGAMEII.*KM.*_a$','PARKGAMEII.*PD.*_a$','PARKGAMEII.*CK.*_c$','PARKGAMEII.*BF.*_c$','PARKGAMEII.*SB.*_a$','PARKGAMEII.*HJ.*_c$'}; %,'PARKGAMEII.*LM.*_a'};
% patient_regex = {'PARKGAMEII.*009_HJ','PARKGAMEII.*013_RP','PARKGAMEII.*027_OR','PARKGAMEII.*046_HJ','PARKGAMEII.*053_LM}; %exclu
for ip = 1 : length(patient_regex)
    clear esuj
    esuj = exam(main_dir,patient_regex{ip});
    esuj.addSerie('ACTIVATION$','run_ACTIVATION',1);
    esuj.getSerie('run_ACTIVATION').addStim(stim_dir, 'MRI_run\d{2}_SPM.mat', 'run', 2 )
    if length(esuj) == 2
       dirFonc(ip,:) = esuj.getSerie('run_ACTIVATION') .toJob;
       stim_files(ip,:) = esuj.getStim.toJob(0);
    end
end

%% old version
%load e % why doesn't the e-object work?!
%%or
e = exam(main_dir,'PARKGAMEII.*_[a,c]$');
%e = exam(main_dir,'PARKGAMEII.*exclu$');
e.addSerie('ACTIVATION$','run_ACTIVATION',1);
e.addSerie('t1mpr.*p2$','anat',1);
[ec,ei] = e.removeIncomplete;
%ei.explore
e = ec;
dir_func_all  = e.getSerie('run_ACTIVATION') .toJob;
dir_anat = e.getSerie('anat').toJob(0); % useful only if individual display in the end (I believe)

%% Make symbolic link of the V2-stim in behav directory of V1-stim

%% Make symbolic link of the V2-wts_OC.nii in ACTIVATION directory of V1-wts
% Make symbolic links from tedana_vtd_mle dir to run dir based on job_meica_afni symbolic link creation 

par.fake = 0;
par.redo = 1;
par.verbose = 2;

%par.subdir        = 'tedana_vtd_mle';
par.subdir        = 'tedana009a1_vtd';
par.sge = 0;
par.run = 1;

% par.warp_file_reg = '^wdn';
% job_symbolic_child_to_parent(dir_func, par);

%par.warp_file_reg = '^s5wts';
%job_symbolic_child_to_parent(dir_func, par);

%par.warp_file_reg = '^s5wdn';
%job_symbolic_child_to_parent(dir_func, par);

par.warp_file_reg = '^s6wts';
job_symbolic_child_to_parent(dir_func_all, par);

%par.warp_file_reg = '^s6wdn';
%job_symbolic_child_to_parent(dir_func, par);

%% make a symbolic link of rp_spm.txt and of multiple_regressors to dir_func
par.subdir = 'wts';
par.regfile_regex = 'multiple_regressors.txt';
for i= 1:length(e)
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
    [job_session(i)] = r_movefile(A_src, A_dst, 'linkn', par);
    job = [job_session];
    
end

%% Create output directories architecture

mkdir(main_dir, model_name{1});
double_model_dir = get_subdir_regex(main_dir, model_name{1});
%double_model_dir = get_subdir_regex(main_dir, 'firstlevel_RS');


%could also use the patient_regex
%patient_list = patient_regex;
patient_list = {'PARKGAMEII_001_NB_a','PARKGAMEII_002_BM_a','PARKGAMEII_003_SM_c','PARKGAMEII_007_SD_a','PARKGAMEII_008_JR_a','PARKGAMEII_023_LJ_c','PARKGAMEII_025_CA_a','PARKGAMEII_028_PC_c','PARKGAMEII_033_DD','PARKGAMEII_039_KM_a','PARKGAMEII_043_PD_a','PARKGAMEII_044_CK_c','PARKGAMEII_047_BF_c','PARKGAMEII_048_SB_a','PARKGAMEII_052_HJ_a','PARKGAMEII_053_LM_c'};
%patient_list = {'P_ARKGAMEII_009_HJ_c','P_ARKGAMEII_013_RP_c','P_ARKGAMEII_027_OR_a','P_ARKGAMEII_046_HJ_c'}; %exclu

for ipatient = 1: length(patient_list)
    mkdir(double_model_dir{1}, patient_list{ipatient});
    patients_dir{ipatient} = get_subdir_regex(double_model_dir, patient_list{ipatient}); %%outdirs
end

%% Job define model
clear par
par.sge = 0;
par.run = 0;
par.display = 0;

par.mask_thr = 0.1;
%par.mask_thr = 0.07;
par.mask = 1; % if we want to use a personalised mask
%par.mask_path = {'/home/anna.skrzatek/data/nifti_test/wmean_mask.nii'}; % specify if we want to use the mean mask, otherwise blank to get the individual warped mask
par.mask_path = {''};
par.TR = 1.6;
par.rp = 1;

%% ts TAPAS

par.file_reg = '^s6wts.*nii';
par.rp_regex = 'wts_multiple_regressors.txt';

par.sge = 1;
par.run = 0;
par.display = 0;
par.jobname = 'spm_glm_auto_dbl_def';

job_auto_doublerun_model_specify(dirFonc, patients_dir, stim_files, par)
%job_ending_rountines(job1,[],par);

%return

%% Model Estimate &

clear par

par.run = 0;
par.sge = 1;
par.display = 0;
par.jobname = 'spm_glm_dble_est'

fspm = addsuffixtofilenames(patients_dir,'SPM.mat');
job_first_level_estimate(fspm,par);

%% Contrast : definition

%% Getting the number of regressors per session - size of the first session contrast //nb of zeros to add

% e.gser('run_ACT').addVolume('wts.*.txt','multireg',1)
regcount1 = [];

for sub = 1 : length(patient_regex)
    clear esuj
    esuj = exam(main_dir,patient_regex{sub});
    esuj.addSerie('ACTIVATION$','run_ACT',1);
    esuj.gser('run_ACT').addVolume('wts.*.txt','multireg',1);
    regpath{sub,:} = esuj.gser('run_ACT').gvol('multireg').path;
%     regpath{sub} = e(sub).gser('run_ACT').gvol('multireg').path;
    regfile = load(regpath{sub,1});
    [nvols nreg] = size(regfile);
    regcount1(sub)= nreg;
    regcon{sub} = zeros(1,nreg);
end
%%

% par.sessrep = 'both';
par.sessrep = 'none';
% par.sessrep = 'replsc';

if par.sessrep == 'none'
    %% T-Contrast : definition
    Rest            = [1 0 0 0 0 0];
    REAL_Left       = [0 1 0 0 0 0];
    REAL_Right      = [0 0 1 0 0 0];
    IMAGINARY_Left  = [0 0 0 1 0 0];
    IMAGINARY_Right = [0 0 0 0 1 0];
    Instruction     = [0 0 0 0 0 1];
    
    for isuj = 1: length(fspm)
        
        %% V1
        Rest_S1            = horzcat(Rest, regcon{isuj},zeros(1,6));
        REAL_Left_S1       = horzcat(REAL_Left, regcon{isuj},zeros(1,6));
        REAL_Right_S1      = horzcat(REAL_Right, regcon{isuj},zeros(1,6));
        IMAGINARY_Left_S1  = horzcat(IMAGINARY_Left, regcon{isuj},zeros(1,6));
        IMAGINARY_Right_S1 = horzcat(IMAGINARY_Right, regcon{isuj},zeros(1,6));
        Instruction_S1     = horzcat(Instruction, regcon{isuj},zeros(1,6));
        
        %% V2
        Rest_S2            = horzcat(zeros(1,6),regcon{isuj},Rest);
        REAL_Left_S2       = horzcat(zeros(1,6),regcon{isuj},REAL_Left);
        REAL_Right_S2      = horzcat(zeros(1,6),regcon{isuj},REAL_Right);
        IMAGINARY_Left_S2  = horzcat(zeros(1,6),regcon{isuj},IMAGINARY_Left);
        IMAGINARY_Right_S2 = horzcat(zeros(1,6),regcon{isuj},IMAGINARY_Right);
        Instruction_S2     = horzcat(zeros(1,6),regcon{isuj},Instruction);
        
        %%
    
        clear contrast_T
            % T contrast

        contrast_T.names = {

            'Rest_S1'
            'Rest_S2'
            'REAL_Left_S1'
            'REAL_Left_S2'
            'REAL_Right_S1'
            'REAL_Right_S2'
            'IMAGINARY_Left_S1'
            'IMAGINARY_Left_S2'
            'IMAGINARY_Right_S1'
            'IMAGINARY_Right_S2'
            'Instruction_S1'
            'Instruction_S2'

            'REAL_Left - Rest_S1'
            'REAL_Left - Rest_S2'
            'REAL_Right - Rest_S1'
            'REAL_Right - Rest_S2'
            'IMAGINARY_Left - Rest_S1'
            'IMAGINARY_Left - Rest_S2'
            'IMAGINARY_Right - Rest_S1'
            'IMAGINARY_Right - Rest_S2'

            %Add the S1 vs S2 contrasts
            'REAL_Right_S2 > S1'
            'IMAGINARY_Right_S2 > S1'
            'REAL_Left_S2 > S1'
            'IMAGINARY_Left_S2 > S1'
            
            'REAL_Right_S1 > S2'
            'IMAGINARY_Right_S1 > S2'
            'REAL_Left_S1 > S2'
            'IMAGINARY_Left_S1 > S2'
            

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

            }';
%%
            contrast_T.values = {

            Rest_S1
            Rest_S2
            REAL_Left_S1
            REAL_Left_S2
            REAL_Right_S1
            REAL_Right_S2
            IMAGINARY_Left_S1
            IMAGINARY_Left_S2
            IMAGINARY_Right_S1
            IMAGINARY_Right_S2
            Instruction_S1
            Instruction_S2

            REAL_Left_S1 - Rest_S1
            REAL_Left_S2 - Rest_S2
            REAL_Right_S1 - Rest_S1
            REAL_Right_S2 - Rest_S2
            IMAGINARY_Left_S1 - Rest_S1
            IMAGINARY_Left_S2 - Rest_S2
            IMAGINARY_Right_S1 - Rest_S1
            IMAGINARY_Right_S2 - Rest_S2

            %Add the S1 vs S2 contrasts
            REAL_Right_S2 - REAL_Right_S1
            IMAGINARY_Right_S2 - IMAGINARY_Right_S1
            REAL_Left_S2 - REAL_Left_S1
            IMAGINARY_Left_S2 - IMAGINARY_Left_S1
            
            REAL_Right_S1 - REAL_Right_S2
            IMAGINARY_Right_S1 - IMAGINARY_Right_S2
            REAL_Left_S1 - REAL_Left_S2
            IMAGINARY_Left_S1 - IMAGINARY_Left_S2

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

            }';


            contrast_T.types = cat(1,repmat({'T'},[1 length(contrast_T.names)]));


            % % F contrast
            % 
            % contrast_F.names = {'F-all'}';
            % 
            % contrast_F.values = {eye(6)}';
            % 
            % 
            % contrast_F.types = cat(1,repmat({'F'},[1 length(contrast_F.names)]));

            % 
            % contrast.names  = [contrast_F.names  contrast_T.names];
            % contrast.values = [contrast_F.values contrast_T.values];
            % contrast.types  = [contrast_F.types  contrast_T.types];

            contrast.names    = [contrast_T.names];
            contrast.values    = [contrast_T.values];
            contrast.types    = [contrast_T.types];

            %% Contrast : write
            clear par

            par.sge = 1;
            par.run = 0;
            par.display = 0;
            par.jobname = sprintf('spm_glm_con_%d',isuj);
            
            par.delete_previous = 1;
            par.report          = 0;

            job_first_level_contrast(fspm(isuj),contrast,par);
    end
end

% %% Display and save results for all subjects, for one contrast example
% 
% e.getSerie('anat').addVolume('^wp0','t1',1);
% t1 = e.getSerie('anat').getVolume('t1');
% Coordlist.values = cat(2,[0.0; 0.0; 0.0], [34.0; -24.0; 68.0], [-34.0; -24.0; 68.0], [-4.0; -48.0; -24.0], [0.0; -24.0; 10.0]); 
% Coordlist.names = {'centre'; 'RIGHT SM'; 'LEFT SM'; 'CEREBELLUM'; 'MOTOR BI'};
% %mkdir('/home/anna.skrzatek/Desktop','auto_figures_first_level_con4');
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

%% Save the new e object

%cd /network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/
%save e

%% Create figures
clear par

%fsub = e.gpath;
fsub = patients_dir;
Coordlist.values = cat(2,[0.0; 0.0; 0.0], [34.0; -24.0; 68.0], [-34.0; -24.0; 68.0], [-4.0; -48.0; -24.0], [0.0; -24.0; 10.0]); % EXAMPLE
Coordlist.names = {'centre'; 'RIGHT SM'; 'LEFT SM'; 'CEREBELLUM'; 'MOTOR BI'}; % EXAMPLE

par.subdir       = 'smodel_tedana';
par.anat_dir_reg = 'S\d{2}_t1mpr_S256_0_8iso_p2$';
par.output_dir   = 'figures';
par.extent       = 5;
par.thresh       = 0.001;
par.fixscale     = 1;
par.minscale     = 0;
par.maxscale     = 20;

%conlist = {'F-all','REAL_Left','REAL_Right','IMAGINARY_Left','IMAGINARY_Right','REAL_Left - Rest','REAL_Right - Rest','IMAGINARY_Left - Rest','IMAGINARY_Right - Rest'};
%conlist = {'REAL_Right - Rest_S1','IMAGINARY_Right - Rest_S1','REAL_Right - Rest_S2','IMAGINARY_Right - Rest_S2'};
conlist = {'REAL_Right_S2 - S1','IMAGINARY_Right_S2 - S1','REAL_Left_S2 - S1','IMAGINARY_Left_S2 - S1'};
addpath /home/anna.skrzatek/matvol/SPM/

for icon = 1 : length(conlist)
   par.conname  = conlist{icon};
   job_spm_single_results_display(fsub, Coordlist, par);
end

close all

% %%script_report_edit(e);

%% Display

%%e.getOne.getModel(model_name{2}).show
