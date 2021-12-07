% Init

clear
clc


addpath /home/anna.skrzatek/matvol/
addpath /home/anna.skrzatek/matvol/SPM/firstlevel/
addpath('/home/anna.skrzatek/MRI_analysis/')

cd      /home/anna.skrzatek/data

main_dir = fullfile(pwd,'/nifti_test');

model_name = {'full_RS_sts_tapas_doublerun_resliced'};%, 'smodel_ts_tapas', 'smodel_dn_tapas'};

%% fetch Input dirs
RSinput_dir = fullfile(main_dir,'/firstlevel_RS')
cd (RSinput_dir)

patient_regex = {'PARKGAMEII.*NB.*_a','PARKGAMEII.*BM.*_a','PARKGAMEII.*SM.*_c','PARKGAMEII.*SD.*_a','PARKGAMEII.*JR.*_a','PARKGAMEII.*LJ.*_c','PARKGAMEII.*CA.*_a','PARKGAMEII.*PC.*_c','PARKGAMEII.*DD.*','PARKGAMEII.*KM.*_a','PARKGAMEII.*PD.*_a','PARKGAMEII.*CK.*_c','PARKGAMEII.*BF.*_c','PARKGAMEII.*SB.*_a','PARKGAMEII.*HJ.*_c'}; %,'PARKGAMEII.*LM.*_c'};
%patient_regex = {'PARKGAMEII_001_NB_a','PARKGAMEII_002_BM_a','PARKGAMEII_003_SM_c','PARKGAMEII_007_SD_a','PARKGAMEII_008_JR_a','PARKGAMEII_023_LJ_c','PARKGAMEII_025_CA_a','PARKGAMEII_028_PC_c','PARKGAMEII_033_DD','PARKGAMEII_039_KM_a','PARKGAMEII_043_PD_a','PARKGAMEII_044_CK_c','PARKGAMEII_047_BF_c','PARKGAMEII_048_SB_a','PARKGAMEII_052_HJ_c'}; 
% patient_regex = {'PARKGAMEII.*009_HJ','PARKGAMEII.*013_RP','PARKGAMEII.*027_OR','PARKGAMEII.*046_HJ','PARKGAMEII.*053_LM}; %exclu

ROIs = {'Caudate_L','Caudate_R','Cereb3_L','Cereb3_R','Cereb4_5_L','Cereb4_5_R','Cereb6_L','Cereb6_R','Cereb7b_L','Cereb7b_R','Cereb8_L','Cereb8_R','Cereb9_L','Cereb9_R','Cereb10_L','Cereb10_R','Cereb_Crus1_L','Cereb_Crus1_R','Cereb_Crus2_L','Cereb_Crus2_R','Cuneus_L','Cuneus_R','Insula_L','Insula_R','Pallidum_L','Pallidum_R','Parietal_Inf_L','Parietal_Inf_R','Parietal_Sup_L','Parietal_Sup_R','PCC','Postcentral_L','Postcentral_R','PPN_L','PPN_R','Precentral_L','Precentral_R','Precuneus_L','Precuneus_R','Putamen_L','Putamen_R','SMA_Ant_L','SMA_Ant_R','SMA_Post_L','SMA_Post_R','Thalamus_L','Thalamus_R','Vermis1_2','Vermis3','Vermis4_5','Vermis6','Vermis7','Vermis8','Vermis9','Vermis10'};
%    ROIs = {'SMA_face_L','SMA_face_R','SMA_foot_L','SMA_foot_R','SMA_hand_L','SMA_hand_R'};
    
dirFonc = cell(length(patient_regex),length(ROIs));

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

        A_src1 = fullfile(reg_src_dir1, par.regfile_regex);
        A_dst1 = fullfile(dirFonc{ip,ir,1}, regfile_out);
        A_src2 = fullfile(reg_src_dir2, par.regfile_regex);
        A_dst2 = fullfile(dirFonc{ip,ir,2}, regfile_out);
    
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

        A_src1 = fullfile(reg_src_dir1, par.regfile_regex);
        A_dst1 = fullfile(dirFonc{ip,ir,1}, par.regfile_regex);
        A_src2 = fullfile(reg_src_dir2, par.regfile_regex);
        A_dst2 = fullfile(dirFonc{ip,ir,2}, regfile_out);
    
            
        par.redo = 0;
        par.verbose = 2;
        par.run = 1;
        par.jobname = sprintf('job_symbolic_link');
        [job_session(i)] = r_movefile(A_src1, A_dst1, 'linkn', par);
        job = [job_session];

        [job_session(i)] = r_movefile(A_src2, A_dst2, 'linkn', par);
        job = [job_session];

    end
end



%% Create output directories architecture

@ -113,11 +99,19 @@ double_model_dir = get_subdir_regex(main_dir, model_name{1});
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



%% Job define model
clear par
par.sge = 0;
par.rp = 1;

%% ts TAPAS

par.file_reg = '^con_0001.nii';
par.rp_regex = 'wts_multiple_regressors.txt';

par.sge = 1;
par.run = 0;
par.display = 0;
par.jobname = 'spm_glm_auto_dbl_def_RS';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TO CREATE a MODIFIED VERSION
job_auto_RS_doublerun_model_specify(dirFonc, outdirs, par)
%job_ending_rountines(job1,[],par);

%return
