%% Analysis for RS groups/ROIs

clc
clear all

main_dir = fullfile('/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek','/nifti_test');
stat_dir = fullfile(main_dir,'full_secondlevel_RS');
A_dir = fullfile(main_dir,'resliced_RS_ANOVA_V1V2');

cd (main_dir)

Stat = exam(stat_dir,'PARK');
%Stat = exam(stat_dir,'PARK_V1');

InputRS   = fullfile(main_dir,'firstlevel_RS'); % newly created - to be checked for transfer errors

patient_list = {'PARKGAMEII_001_NB_a','PARKGAMEII_002_BM_a','PARKGAMEII_003_SM_c','PARKGAMEII_007_SD_a','PARKGAMEII_008_JR_a','PARKGAMEII_023_LJ_c','PARKGAMEII_025_CA_a','PARKGAMEII_028_PC_c','PARKGAMEII_033_DD','PARKGAMEII_039_KM_a','PARKGAMEII_043_PD_a','PARKGAMEII_044_CK_c','PARKGAMEII_047_BF_c','PARKGAMEII_048_SB_a','PARKGAMEII_052_HJ_c'};
patient_list_V1 = {'PARKGAMEII_001_NB_a','PARKGAMEII_002_BM_a','PARKGAMEII_003_SM_c','PARKGAMEII_007_SD_a','PARKGAMEII_008_JR_a','PARKGAMEII_023_LJ_c','PARKGAMEII_025_CA_a','PARKGAMEII_028_PC_c','PARKGAMEII_033_DD','PARKGAMEII_039_KM_a','PARKGAMEII_040_RE_a','PARKGAMEII_042_RS_a','PARKGAMEII_043_PD_a','PARKGAMEII_044_CK_c','PARKGAMEII_046_HJ_c','PARKGAMEII_047_BF_c','PARKGAMEII_048_SB_a','PARKGAMEII_052_HJ_c','PARKGAMEII_053_LM_a'};

    for ipatient = 1: length(patient_list)
        mkdir(InputRS, patient_list{ipatient});
        patients_dir{ipatient} = get_subdir_regex(InputRS, patient_list{ipatient}); %%outdirs
        patients_dir_V1{ipatient} = get_subdir_regex(InputRS, patient_list_V1{ipatient}); %%outdirs
    end
    
RSObj_a = exam(InputRS,'PARK.*a$');
RSObj_c = exam(InputRS,'PARK.*[c,DD]$');
RS_all = exam( InputRS,'PARK.*a$') + exam(InputRS,'PARK.*[c,DD]$');
RS_V1_all = exam( InputRS,'PARK.*a$') + exam(InputRS,'PARK.*a_V1') + exam(InputRS,'PARK.*[c,DD]$') + exam(InputRS,'PARK.*c_V1');

ROIs = {'Caudate_L','Caudate_R','Cereb3_L','Cereb3_R','Cereb4_5_L','Cereb4_5_R','Cereb6_L','Cereb6_R','Cereb7b_L','Cereb7b_R','Cereb8_L','Cereb8_R','Cereb9_L','Cereb9_R','Cereb10_L','Cereb10_R','Cereb_Crus1_L','Cereb_Crus1_R','Cereb_Crus2_L','Cereb_Crus2_R','Cuneus_L','Cuneus_R','Insula_L','Insula_R','Pallidum_L','Pallidum_R','Parietal_Inf_L','Parietal_Inf_R','Parietal_Sup_L','Parietal_Sup_R','PCC','Postcentral_L','Postcentral_R','PPN_L','PPN_R','Precentral_L','Precentral_R','Precuneus_L','Precuneus_R','Putamen_L','Putamen_R','SMA_Ant_L','SMA_Ant_R','SMA_Face_L','SMA_Face_R','SMA_Foot_L','SMA_Foot_R','SMA_Hand_L','SMA_Hand_R','SMA_Post_L','SMA_Post_R','Thalamus_L','Thalamus_R','Vermis1_2','Vermis3','Vermis4_5','Vermis6','Vermis7','Vermis8','Vermis9','Vermis10'};
%ROIs = {'SMA_Ant_L','SMA_Ant_R','SMA_Face_L','SMA_Face_R','SMA_Foot_L','SMA_Foot_R','SMA_Hand_L','SMA_Hand_R','SMA_Post_L','SMA_Post_R'};
%ROIs = {'Cereb3_L','Cereb3_R','Cereb4_5_L','Cereb4_5_R','Cereb6_L','Cereb6_R','Cereb7b_L','Cereb7b_R','Cereb8_L','Cereb8_R','Cereb9_L','Cereb9_R','Cereb10_L','Cereb10_R','Cereb_Crus1_L','Cereb_Crus1_R','Cereb_Crus2_L','Cereb_Crus2_R'};
Conditions  = ROIs;
Sessions    = {'V1','V2'}%,'V2-V1'};
A_Sessions  = {'V2_V1'}; % pour ANOVA

for iR = 1 : length(ROIs)
    for iS = 1 : length(Sessions)
        %mkdir(StatDir{:},char(sprintf('%s_%s',Conditions{iC},Sessions{iS})))
        RSObj_a.mkdir(sprintf('%s_%s',ROIs{iR},Sessions{iS}));
        RSObj_a.addSerie(sprintf('%s_%s',ROIs{iR},Sessions{iS}),sprintf('%s_%s',ROIs{iR},Sessions{iS}));

        % GROUP A add volumes for each contrast - each individual
        RSObj_a.getSerie(sprintf('%s_%s',ROIs{iR},Sessions{iS})).addVolume('con.*01',sprintf('%s_%s',ROIs{iR},Sessions{iS}))

        RSObj_c.mkdir(sprintf('%s_%s',ROIs{iR},Sessions{iS}));
        RSObj_c.addSerie(sprintf('%s_%s',ROIs{iR},Sessions{iS}),sprintf('%s_%s',ROIs{iR},Sessions{iS}));

        % GROUP C add volumes for each contrast - each individual
        RSObj_c.getSerie(sprintf('%s_%s',ROIs{iR},Sessions{iS})).addVolume('con.*01',sprintf('%s_%s',ROIs{iR},Sessions{iS}))
        
    end
    mkdir(A_dir,sprintf('%s_V2_V1',ROIs{iR}))
    
    % FIRSTLEVEL for ALL PARKs add volumes for each ROI's contrast
     
     RS_V1_all.addSerie(sprintf('%s_V1$',ROIs{iR}),sprintf('%s_V1',ROIs{iR}),1);
     RS_V1_all.getSerie(sprintf('%s_V1$',ROIs{iR})).addVolume('con.*01',sprintf('%s_V1_con$',ROIs{iR}))
end


clear RS_allcell gRS_all

n = 1;
for iC = 1 : length(Conditions)
    for iS = 1 : length(Sessions)
        %mkdir(StatDir{:},char(sprintf('%s_%s',Conditions{iC},Sessions{iS})))
        Stat.mkdir(sprintf('%s_%s',Conditions{iC},Sessions{iS}));
        Stat.addSerie(sprintf('^%s_%s$',Conditions{iC},Sessions{iS}),sprintf('%s_%s',Conditions{iC},Sessions{iS}));
        
        % cons a & cons_c to be found with the
        % sprintf('%s_%s',Conditions{iC},Sessions{iS}) match and a
        % structure needed for the output of it (not to lose any cons by replacing them)
        % we can use the getSerie.getVol .toJob with the regex of the above 

        RS_a{n}   = RSObj_a.getSerie(sprintf('%s_%s',Conditions{iC},Sessions{iS})).getVolume(sprintf('%s_%s',Conditions{iC},Sessions{iS})) %.toJob;
        gRS_a{n} = cellstr(RS_a{n}(:) .toJob)            
        RS_c{n}  = RSObj_c.getSerie(sprintf('%s_%s',Conditions{iC},Sessions{iS})).getVolume(sprintf('%s_%s',Conditions{iC},Sessions{iS})) %.toJob;
        gRS_c{n} = cellstr(RS_c{n}(:) .toJob)
        
        n = n + 1;
    end
    RS_allcell{iC} = RS_V1_all.getSerie(sprintf('%s_V1',Conditions{iC})).getVolume(sprintf('%s_V1',Conditions{iC})) %.toJob;
    gRS_all{iC} = cellstr(RS_allcell{iC}(:) .toJob)        
end


%outdirs = Stat.getSerie('Cereb.*V1') .toJob; % 7 x 1 x 108 cells
outdirs = Stat.getSerie('.*V') .toJob; % 7 x 1 x 108 cells
outdirs = outdirs(2:3);
A_outdir = exam(A_dir,'V2_V1').toJob;
groups  = {gRS_a, gRS_c};
outdirs_V1 = Stat(1).getSerie('.*V1') .toJob; % 7 x 1 x 108 cells
group_V1 = gRS_all;

%% Regressors definition
%% AGE V2
covars{1,1} = [70 
74 
64 
76 
79 
61 
75 
66];

covars{1,2} = [72 
68 
68 
72 
56 
57 
66];

%% GENDER V2 % F : 1 M : 2
covars{2,1} = [1
1
1
2
1
1
2
2];

covars{2,2} = [1
2
2
2
2
1
2];

%% AGE V1
covars_V1{1,1} = [70 
74 
64 
76 
79 
61 
75 
66
71
72
59];

covars_V1{1,2} = [72 
68 
68 
72 
56 
57 
66 
62];

%% GENDER V1 % F : 1 M : 2
covars_V1{2,1} = [1
1
1
2
1
1
2
2
1
2
1];

covars_V1{2,2} = [1
2
2
2
2
1
2
1];

%% per group
addpath /home/anna.skrzatek/MRI_analysis/

par.run = 0;
par.sge = 1;
par.jobname = 'job_RS_secondlevel_auto_ac';

secondlevel_RS_matlabbatch(groups,outdirs,covars,par)

%% V1 session only

par.run = 1;
par.sge = 0;
par.jobname = 'job_RS_V1_secondlevel_V1_auto';

all_covars{1,1} = vertcat(covars_V1{1,1},covars_V1{1,2});
all_covars{2,1} = vertcat(covars_V1{2,1},covars_V1{2,2});
groups_combined = {group_V1};

secondlevel_RS_matlabbatch(groups_combined,outdirs_V1,all_covars,par)

%% models estimate

clear par
cd (main_dir)

for iout = 1 : length(outdirs)
    
    fspm = addsuffixtofilenames(outdirs{iout}, 'SPM.mat');

    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
    par.jobname  = 'spm_secondlevel_RS_est';
    job_first_level_estimate(fspm,par)
    
end

%% V1 only

clear par
cd (main_dir)

for iout = 1 : length(outdirs_V1)
    
    fspm = addsuffixtofilenames(outdirs_V1{iout}, 'SPM.mat');

    par.run = 1;
    %par.sge = 1;
    par.sge_queu = 'normal,bigmem';
    par.jobname  = 'spm_secondlevel_RS_est';
    job_first_level_estimate(fspm,par)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Contrast creation for each SPM.mat

% t-statistics
    Main = 1;

%%

for iout = 1 : length(outdirs)
    %fspm = addsuffixtofilenames(outdirs{iout},'SPM.mat');
    modest = addsuffixtofilenames(outdirs{iout},'SPM.mat');
    for iroi = 1 : length(modest)
        parts = strsplit(char(modest(iroi)), '/');
        roilabel = parts{end-1};   
    
        %% Contrast names
        contrast_t.names = {
            sprintf('Connectivity %s',roilabel)
            }';

        %% Contrast values
        contrast_t.values = {
            Main
            }';

        %% Contrast type
        contrast_t.types = cat(1,repmat({'T'},[1 length(contrast_t.names)]));

        contrast.names  = [contrast_t.names];
        contrast.values = [contrast_t.values];
        contrast.types  = [contrast_t.types];

        %% Contrast : write
        clear par

        par.sge = 0;
        par.run = 1;
        par.display = 0;
        par.jobname = sprintf('spm_write_%s_con',roilabel);

        % par.sessrep = 'both';
        par.sessrep = 'none';

        par.delete_previous = 1;
        par.report          = 0;

        job_first_level_contrast(modest(iroi),contrast,par);
        
%         Stat(iout).getSerie(roilabel).addVolume('spmT_0001','main',1)
        
%         mainef = Stat(iout).getSerie(roilabel).getVolume('main') .toJob
        
        mask{iroi} = cellstr(fullfile(outdirs{iout}{iroi},'mask.nii'));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% V1
% t-statistics
    Main = 1;

%%

for iout = 1 : length(outdirs_V1)
    %fspm = addsuffixtofilenames(outdirs{iout},'SPM.mat');
    modest = addsuffixtofilenames(outdirs_V1{iout},'SPM.mat');
    for iroi = 1 : length(modest)
        parts = strsplit(char(modest(iroi)), '/');
        roilabel = parts{end-1};   
    
        %% Contrast names
        contrast_t.names = {
            sprintf('Connectivity %s',roilabel)
            }';

        %% Contrast values
        contrast_t.values = {
            Main
            }';

        %% Contrast type
        contrast_t.types = cat(1,repmat({'T'},[1 length(contrast_t.names)]));

        contrast.names  = [contrast_t.names];
        contrast.values = [contrast_t.values];
        contrast.types  = [contrast_t.types];

        %% Contrast : write
        clear par

        par.sge = 0;
        par.run = 1;
        par.display = 0;
        par.jobname = sprintf('spm_write_%s_con_V1',roilabel);

        % par.sessrep = 'both';
        par.sessrep = 'none';

        par.delete_previous = 1;
        par.report          = 0;

        job_first_level_contrast(modest(iroi),contrast,par);
        
%         Stat(iout).getSerie(roilabel).addVolume('spmT_0001','main',1)
        
%         mainef = Stat(iout).getSerie(roilabel).getVolume('main') .toJob
        
        mask{iroi} = cellstr(fullfile(outdirs_V1{iout}{iroi},'mask.nii'));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% paired t-test models definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd (main_dir)

% define outdirs

StatDp = exam(stat_dir,'deltas','.*paired');
Sessions    = {'V2_V1'};%

n = 1;
for iC = 1 : length(Conditions)
    for iS = 1 : length(Sessions)
        StatDp.mkdir(sprintf('%s_%s',Conditions{iC},Sessions{iS}));
        StatDp.addSerie(sprintf('^%s_%s$',Conditions{iC},Sessions{iS}),sprintf('%s_%s',Conditions{iC},Sessions{iS}));
        n = n + 1;
        
    end
end

outdirs = StatDp.getSerie('.*') .toJob; % 2 x 55 x 1 cells

% define the volumes per group
groups  = {gRS_a, gRS_c};

for igroup = 1 : length(groups)
    i = 0;
    for icov = 1: length(covars{1,igroup})
        i = i+1;
        covars_double{1,igroup}(i) = covars{1,igroup}(icov);
        covars_double{2,igroup}(i) = covars{2,igroup}(icov);
        i = i+1;
        covars_double{1,igroup}(i) = covars{1,igroup}(icov);
        covars_double{2,igroup}(i) = covars{2,igroup}(icov);    
    end
end

% define models per group
addpath /home/anna.skrzatek/MRI_analysis/

par.run = 0;
par.sge = 1;
par.jobname = 'job_RS_paired_secondlevel_auto';

secondlevel_paired_RS_matlabbatch(groups,outdirs,covars_double,par)

%% models estimate
clear par
cd (main_dir)

for iout = 1 : length(outdirs)
    
    fspm = addsuffixtofilenames(outdirs{iout}, 'SPM.mat');

    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
    par.jobname  = 'spm_paired_secondlevel_RS_est';
    job_first_level_estimate(fspm,par)
    
end

%% Contrast creation for each SPM.mat

% t-statistics
    V1_V2 = [ 1 -1];
    V2_V1 = [-1  1];

for iout = 1 : length(outdirs)
    %fspm = addsuffixtofilenames(outdirs{iout},'SPM.mat');
    modest = addsuffixtofilenames(outdirs{iout},'SPM.mat');
    for iroi = 1 : length(modest)
        parts = strsplit(char(modest(iroi)), '/');
        roilabel = parts{end-1};   
    
        %% Contrast names
        contrast_t.names = {
            sprintf('Connectivity V1>V2 %s',roilabel)
            sprintf('Connectivity V1<V2 %s',roilabel)
            }';

        %% Contrast values
        contrast_t.values = {
            V1_V2
            V2_V1
            }';

        %% Contrast type
        contrast_t.types = cat(1,repmat({'T'},[1 length(contrast_t.names)]));

        contrast.names  = [contrast_t.names];
        contrast.values = [contrast_t.values];
        contrast.types  = [contrast_t.types];

        %% Contrast : write
        clear par

        par.sge = 0;
        par.run = 1;
        par.display = 0;
        par.jobname = sprintf('spm_write_%s_con',roilabel);

        % par.sessrep = 'both';
        par.sessrep = 'none';

        par.delete_previous = 1;
        par.report          = 0;

        job_first_level_contrast(modest(iroi),contrast,par);
        
%         Stat(iout).getSerie(roilabel).addVolume('spmT_0001','main',1)
        
%         mainef = Stat(iout).getSerie(roilabel).getVolume('main') .toJob
        
        mask{iroi} = cellstr(fullfile(outdirs{iout}{iroi},'mask.nii'));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2 sample t-test models definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd (main_dir)

% define outdirs

StatD2s = exam(stat_dir,'deltas','.*2sample');
Sessions    = {'V1','V2'};%

n = 1;
for iC = 1 : length(Conditions)
    for iS = 1 : length(Sessions)
        StatD2s.mkdir(sprintf('%s_%s',Conditions{iC},Sessions{iS}));
        StatD2s.addSerie(sprintf('^%s_%s$',Conditions{iC},Sessions{iS}),sprintf('%s_%s',Conditions{iC},Sessions{iS}));
        n = n + 1;
        
    end
end

outdirs = StatD2s.getSerie('.*') .toJob; % 55 x 1 cells

% define the volumes per group
groups  = {gRS_a, gRS_c};

% define models per group
addpath /home/anna.skrzatek/MRI_analysis/

par.run = 0;
par.sge = 1;
par.jobname = 'job_RS_2s_secondlevel_auto_covars';

secondlevel_RS_2s_matlabbatch(groups,outdirs,covars,par)

%% models estimate
clear par
cd (main_dir)

for iout = 1 : length(outdirs)
    
    fspm = addsuffixtofilenames(outdirs{iout}, 'SPM.mat');

    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
    par.jobname  = 'spm_secondlevel_RS_2s_est';
    job_first_level_estimate(fspm,par)
    
end

%% Contrast creation for each SPM.mat

% t-statistics
    Kinect_Ordi = [ 1 -1];
    Ordi_Kinect = [-1  1];
    
for iout = 1 : length(outdirs)
    %fspm = addsuffixtofilenames(outdirs{iout},'SPM.mat');
    modest = addsuffixtofilenames(outdirs{iout},'SPM.mat');
    for iroi = 1 : length(modest)
        parts = strsplit(char(modest(iroi)), '/');
        roilabel = parts{end-1};   
    
        %% Contrast names
        contrast_t.names = {
            sprintf('Connectivity K > O %s',roilabel)
            sprintf('Connectivity K < O %s',roilabel)
            }';

        %% Contrast values
        contrast_t.values = {
            Kinect_Ordi
            Ordi_Kinect
            }';

        %% Contrast type
        contrast_t.types = cat(1,repmat({'T'},[1 length(contrast_t.names)]));

        contrast.names  = [contrast_t.names];
        contrast.values = [contrast_t.values];
        contrast.types  = [contrast_t.types];

        %% Contrast : write
        clear par

        par.sge = 0;
        par.run = 1;
        par.display = 0;
        par.jobname = sprintf('spm_write_%s_con',roilabel);

        % par.sessrep = 'both';
        par.sessrep = 'none';

        par.delete_previous = 1;
        par.report          = 0;

        job_first_level_contrast(modest(iroi),contrast,par);
        
%         Stat(iout).getSerie(roilabel).addVolume('spmT_0001','main',1)
        
%         mainef = Stat(iout).getSerie(roilabel).getVolume('main') .toJob
        
        mask{iroi} = cellstr(fullfile(outdirs{iout}{iroi},'mask.nii'));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% full factorial ANOVA model definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 1;
for iC = 1 : length(Conditions)
    RS_a_V1{n}   = RSObj_a.getSerie(sprintf('%s_V1',Conditions{iC})).getVolume(sprintf('%s_V1',Conditions{iC})) %.toJob;
    gRS_a_V1{n} = cellstr(RS_a_V1{n}(:) .toJob)            
    RS_c_V1{n}  = RSObj_c.getSerie(sprintf('%s_V1',Conditions{iC})).getVolume(sprintf('%s_V1',Conditions{iC})) %.toJob;
    gRS_c_V1{n} = cellstr(RS_c_V1{n}(:) .toJob)

    RS_a_V2{n}   = RSObj_a.getSerie(sprintf('%s_V2',Conditions{iC})).getVolume(sprintf('%s_V2',Conditions{iC})) %.toJob;
    gRS_a_V2{n} = cellstr(RS_a_V2{n}(:) .toJob)            
    RS_c_V2{n}  = RSObj_c.getSerie(sprintf('%s_V2',Conditions{iC})).getVolume(sprintf('%s_V2',Conditions{iC})) %.toJob;
    gRS_c_V2{n} = cellstr(RS_c_V2{n}(:) .toJob)
    
    n = n+1;
end

% define the volumes per group
A_groups  = {gRS_a_V1, gRS_a_V2; gRS_c_V1, gRS_c_V2}; % A_groups{1,:} : 2 groupes en V1 ; A_groups{:,1} : 2 visites de Kinect

% define models
cd (main_dir)
addpath /home/anna.skrzatek/MRI_analysis/

par.run = 0;
par.sge = 1;
par.jobname = 'job_RS_ANOVA_full_auto_covars';
par.con_auto = 0;

ANOVA_RS_matlabbatch(A_groups,A_outdir,covars,par)

%% models estimate
clear par
cd (main_dir)

for iout = 1 : length(A_outdir)
    
    fspm = addsuffixtofilenames(A_outdir(iout), 'SPM.mat');

    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
    par.jobname  = 'spm_ANOVA_RS_est';
    job_first_level_estimate(fspm,par)
    
end

%% Contrast creation for each SPM.mat

% F-statistics
    Average_condition_effect         = [ 1 1  1  1 0 0];
    Main_Session_effect              = [ 1 1 -1 -1 0 0];
    Main_Group_effect                = [ 1 -1 1 -1 0 0];
    Main_Interaction_SxG_effect      = [ 1 -1 -1 1 0 0];
    
% t-statistics
    Positive_condition_eff           = [1 1 1 1 0 0];
    Negative_condition_eff           = [-1 -1 -1 -1 0 0];
    Positive_Training_eff            = [-1 -1 1 1 0 0];
    Negative_Training_eff            = [1 1 -1 -1 0 0];
    Positive_Kinect_eff              = [1 -1 1 -1 0 0];
    Positive_Ordi_eff                = [-1 1 -1 1 0 0];
    Positive_Interaction_eff         = [1 -1 -1 1 0 0];
    Negative_Interaction_eff         = [-1 1 1 -1 0 0];

%% for all ROIs

for iout = 1 : length(A_outdir)
    %fspm = addsuffixtofilenames(A_outdir{iout},'SPM.mat');
    modest = addsuffixtofilenames(A_outdir{iout},'SPM.mat');
    parts = strsplit(modest, '/');
    roilabel = parts{end-1};   

    %% Contrast names
    contrast_F.names = {
        sprintf('Average condition effect %s',roilabel)
        sprintf('Main session effect %s',roilabel)
        sprintf('Main group effect %s',roilabel)
        sprintf('Main SxG Interaction effect %s',roilabel)
        };

    contrast_t.names = {
        sprintf('Positive Condition effect %s',roilabel)
        sprintf('Negative Condition effect %s',roilabel)
        sprintf('Positive Training effect %s',roilabel)
        sprintf('Negative Training effect %s',roilabel)
        sprintf('Positive Kinect effect %s',roilabel)
        sprintf('Positive Ordi effect %s',roilabel)
        sprintf('Positive Interaction effect %s',roilabel)
        sprintf('Negative Interaction effect %s',roilabel)
        }';

    %% Contrast values
    contrast_F.values = {
        Average_condition_effect
        Main_Session_effect
        Main_Group_effect
        Main_Interaction_SxG_effect
        }';

    contrast_t.values = {
        Positive_condition_eff
        Negative_condition_eff
        Positive_Training_eff
        Negative_Training_eff
        Positive_Kinect_eff
        Positive_Ordi_eff
        Positive_Interaction_eff
        Negative_Interaction_eff
        }';

    %% Contrast type
    contrast_F.types = cat(1,repmat({'F'},[1 length(contrast_F.names)]));
    contrast_t.types = cat(1,repmat({'T'},[1 length(contrast_t.names)]));

%    contrast.names  = [contrast_F.names contrast_t.names];
%    contrast.values = [contrast_F.values contrast_t.values];
%    contrast.types  = [contrast_F.types contrast_t.types];

    %% F-Contrast : write
    clear par

    par.sge = 0;
    par.run = 1;
    par.display = 0;
    par.jobname = sprintf('spm_write_%s_F_con',roilabel);

    % par.sessrep = 'both';
    par.sessrep = 'none';

    par.delete_previous = 1;
    par.report          = 0;

    contrast.names  = [contrast_F.names];
    contrast.values = [contrast_F.values];
    contrast.types  = [contrast_F.types];

    job_first_level_contrast({modest},contrast,par);

    %% t-Contrast : write

    par.jobname = sprintf('spm_write_%s_t_con',roilabel);
    par.delete_previous = 0;
    
    contrast.names  = [contrast_t.names];
    contrast.values = [contrast_t.values];
    contrast.types  = [contrast_t.types];
    
    job_first_level_contrast({modest},contrast,par);

%         Stat(iout).getSerie(roilabel).addVolume('spmT_0001','main',1)

%         mainef = Stat(iout).getSerie(roilabel).getVolume('main') .toJob

%    mask{iroi} = cellstr(fullfile(A_outdir{iout}{iroi},'mask.nii'));
end


