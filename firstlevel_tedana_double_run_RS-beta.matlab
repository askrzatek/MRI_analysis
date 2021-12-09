%% adapted from script RS firstlevel : Cecile Gallea 
%% Init

clear
clc

addpath /home/anna.skrzatek/matvol/SPM/firstlevel/
addpath('/home/anna.skrzatek/MRI_analysis/')

main_dir = '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test';
cd(main_dir)

%addpath '/network/lustre/iss01/cenir/analyse/irm/users/cecile.gallea/ASYA/asyasuit/'

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
mkdir(model_name);

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

i = 1;
for ip = 1 : length(patient_regex)    
    clear esuj
    esuj = exam(RSinput_dir,patient_regex{ip});
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
        regfile_out = sprintf('%s_%s',par.subdir,par.regfile_regex);

        A_src1 = fullfile(reg_src_dir1{:}, par.regfile_regex);
        A_src2 = fullfile(reg_src_dir2{:}, par.regfile_regex);
        A_dst = fullfile(dirFonc{ip,ir}, regfile_out);
        A_dst1 = A_dst{1};
        A_dst2 = A_dst{2};

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

        A_src1 = fullfile(reg_src_dir1{:}, par.regfile_regex);
        A_src2 = fullfile(reg_src_dir2{:}, par.regfile_regex);
        A_dst = fullfile(dirFonc{ip,ir}, regfile_out);
        A_dst1 = A_dst{1};
        A_dst2 = A_dst{2};

            
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
        regfile_out1 = regfile_out1(end-11:end);
        A_src2 = get_subdir_regex_files(reg_src_dir2{:}, par.regfile_regex);
        regfile_out2 = char(A_src1);
        regfile_out2 = regfile_out2(end-11:end);
        A_dst = fullfile(dirFonc{ip,ir}, regfile_out);
        A_dst1 = A_dst{1};
        A_dst2 = A_dst{2};

        par.redo = 0;
        par.verbose = 2;
        par.run = 1;
        par.jobname = sprintf('job_symbolic_link');
        [job_session(i)] = r_movefile(A_src1, A_dst1, 'linkn', par);
        job = [job_session];

        [job_session(i)] = r_movefile(A_src2, A_dst2, 'linkn', par);
        job = [job_session];

    end

    e.getSerie('RS').addVolume('^s6wts','s6wts',1);
    e.getSerie('RS').addVolume('wts_multiple_regressors','wts_rp',1);
    e.getSerie('RS').addVolume('^wbet.*mask','wbet_mask',1);

    RS_all = RS_all + e;

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
end

clear y1 y2 y3 y4 y5 y6 y7 y8 t Matrix Mat_sub u1 u2 u3 u4 u5 u6 u7 u8



%% Parameters
clear par

par.TR = 1.6;

par.file_reg = '^s6wts.*nii';
par.rp       = 1;
par.rp_regex = 'wts_multiple.*txt';

% Masking
par.mask_thr = 0.1; % spm default option
%par.mask     =  {}; % cell(char) of the path for the mask of EACH model : N models means N paths
par.mask     =  gpath(e.gser('run_RS').gvol('wmask')); % cell(char) of the path for the mask of EACH model : N models means N paths
par.cvi      = 'AR(1)'; % 'AR(1)' / 'FAST' / 'none'

% Regressors
%-----------
% multilevel_cells(struct) for used defined regressors : they will NOT be convolved
par.user_regressor = {};
% multilevel_cells(char  ) for used defined regressors : they will NOT be convolved
% The regressors in the file will be concatenated with rp_*.txt
par.file_regressor = addsuffixtofilenames(dirFonc, '/wts_multiple_regressors.txt');

par.jobname  = 'spm_glm_rs_wbet';
par.walltime = '04:00:00';

par.sge   = 0;
par.run      = 0;
par.display  = 1;
par.redo     = 0;





%% Job define model
clear par

par.sge = 0;
par.run = 1;
par.display = 0;

par.mask_thr = 0.1;
par.mask = 1; % if we want to use a personalised mask
par.mask_path = {''};
par.TR = 1.6;
par.rp = 1;
par.file_reg = '^con_0001.nii';
par.rp_regex = 'wts_multiple_regressors.txt';

par.jobname = 'spm_glm_auto_dbl_def_RS';

job_RS_doublerun_model_specify(dirFonc, outdirs, par)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    esuj.addSerie('RS$','run_ACT',1);
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
            'Instruction_S1'

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create output directories architecture
patient_list = {'PARKGAMEII_001_NB_a','PARKGAMEII_002_BM_a','PARKGAMEII_003_SM_c','PARKGAMEII_007_SD_a','PARKGAMEII_008_JR_a','PARKGAMEII_023_LJ_c','PARKGAMEII_025_CA_a','PARKGAMEII_028_PC_c','PARKGAMEII_033_DD','PARKGAMEII_039_KM_a','PARKGAMEII_043_PD_a','PARKGAMEII_044_CK_c','PARKGAMEII_047_BF_c','PARKGAMEII_048_SB_a','PARKGAMEII_052_HJ_c'}; 
%patient_list = {'P_ARKGAMEII_009_HJ_c','P_ARKGAMEII_013_RP_c','P_ARKGAMEII_027_OR_a','P_ARKGAMEII_046_HJ_c','PARKGAMEII_053_LM_c'}; %exclu

Conditions = ROIs;

for ipatient = 1: length(patient_list)
    mkdir(double_model_dir{1}, patient_list{ipatient});
    patients_dir{ipatient} = get_subdir_regex(double_model_dir, patient_list{ipatient});
    for iROI = 1: length(Conditions)
        mkdir(char(patients_dir{ipatient}),Conditions{iROI});
        outdirs{ipatient,iROI} =  get_subdir_regex(patients_dir{ipatient},Conditions{iROI});%% outdirs
    end
end

close all

% %%script_report_edit(e);

%% Display

%%e.getOne.getModel(model_name{2}).show
