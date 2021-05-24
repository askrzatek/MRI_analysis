%% Analysis for RS groups/ROIs

clc
clear all

main_dir = fullfile('/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek','/nifti_test');
stat_dir = fullfile(main_dir,'resliced_secondlevel_RS');

cd (main_dir)

Stat = exam(stat_dir,'PARK');

InputRS   = fullfile(main_dir,'firstlevel_RS'); % newly created - PCC models missing

patient_list = {'PARKGAMEII_001_NB_a','PARKGAMEII_002_BM_a','PARKGAMEII_003_SM_c','PARKGAMEII_007_SD_a','PARKGAMEII_008_JR_a','PARKGAMEII_023_LJ_c','PARKGAMEII_025_CA_a','PARKGAMEII_028_PC_c','PARKGAMEII_033_DD','PARKGAMEII_039_KM_a','PARKGAMEII_043_PD_a','PARKGAMEII_044_CK_c','PARKGAMEII_047_BF_c','PARKGAMEII_048_SB_a'};

    for ipatient = 1: length(patient_list)
        mkdir(InputRS, patient_list{ipatient});
        patients_dir{ipatient} = get_subdir_regex(InputRS, patient_list{ipatient}); %%outdirs
    end
    
RSObj_a = exam(InputRS,'PARK.*a$');
RSObj_c = exam(InputRS,'PARK.*[c,DD]$');

ROIs = {'Caudate_L','Caudate_R','Cereb3_L','Cereb3_R','Cereb4_5_L','Cereb4_5_R','Cereb6_L','Cereb6_R','Cereb7b_L','Cereb7b_R','Cereb8_L','Cereb8_R','Cereb9_L','Cereb9_R','Cereb10_L','Cereb10_R','Cereb_Crus1_L','Cereb_Crus1_R','Cereb_Crus2_L','Cereb_Crus2_R','Cuneus_L','Cuneus_R','Insula_L','Insula_R','Pallidum_L','Pallidum_R','Parietal_Inf_L','Parietal_Inf_R','Parietal_Sup_L','Parietal_Sup_R','PCC','Postcentral_L','Postcentral_R','PPN_L','PPN_R','Precentral_L','Precentral_R','Precuneus_L','Precuneus_R','Putamen_L','Putamen_R','SMA_Ant_L','SMA_Ant_R','SMA_Post_L','SMA_Post_R','Thalamus_L','Thalamus_R','Vermis1_2','Vermis3','Vermis4_5','Vermis6','Vermis7','Vermis8','Vermis9','Vermis10'};
Conditions  = ROIs;
Sessions    = {'V1','V2'}%,'V2-V1'};

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
end

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
end


outdirs = Stat.getSerie('.*') .toJob; % 7 x 1 x 108 cells
groups  = {gRS_a, gRS_c};

%% Regressors definition
%% AGE
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
57];

%% GENDER
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
1];

addpath /home/anna.skrzatek/MRI_analysis/

par.run = 0;
par.sge = 1;
par.jobname = 'job_RS_secondlevel_auto';

secondlevel_RS_matlabbatch(groups,outdirs,covars,par)

%% models estimate

clear par
cd (main_dir)

for iout = 1 : length(outdirs)
    
    fspm = addsuffixtofilenames(outdirs{iout}, 'SPM.mat');

    par.run = 1;
    %par.sge = 1;
    par.sge_queu = 'normal,bigmem';
    par.jobname  = 'spm_secondlevel_RS_est';
    job_first_level_estimate(fspm,par)
    
end

%% Contrast creation for each SPM.mat

% t-statistics
    Main = 1;

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

cd (main_dir)

StatD = exam(stat_dir,'deltas');
Sessions    = {'V2_V1'}%

n = 1;
for iC = 1 : length(Conditions)
    for iS = 1 : length(Sessions)
        StatD.mkdir(sprintf('%s_%s',Conditions{iC},Sessions{iS}));
        StatD.addSerie(sprintf('^%s_%s$',Conditions{iC},Sessions{iS}),sprintf('%s_%s',Conditions{iC},Sessions{iS}));
        n = n + 1;
        
    end
end

outdirs = StatD.getSerie('.*') .toJob; % 7 x 1 x 108 cells
groups  = {gRS_a, gRS_c};

addpath /home/anna.skrzatek/MRI_analysis/

par.run = 0;
par.sge = 1;
par.jobname = 'job_RS_secondlevel_auto';

secondlevel_RS_matlabbatch(groups,outdirs,par)

