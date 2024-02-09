%% paramètres / variables
clc
clear

main_dir = '/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/DATA/Non_chirurgicaux';
project_dir = '/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG';
cd (project_dir)
outdirs = {'/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Results/RS_2sample_ttest/ALFF','/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Results/RS_2sample_ttest/ALFF'} ;

load('e_nonchir.mat');
% fichier de correspondance numero IRM - comportement - groupe - age
d = readtable( [ './DATA/' , 'Correspondance_Numero_Comportement_IRM.csv' ])  ;

%% Regressors definition
% importing the table with multiple columns with patients characteristics
% from CSV: filtered by IRM==1
tab = readtable(fullfile(project_dir,'DATA/ANT_Alerting_RT_multiregression.csv'));
%tab_sel = ismember(tab.IRM, 1);
%group = tab.Group(tab_sel);
%age   = tab.Age(tab_sel);
%sex   = tab.Genre(tab_sel);
%ANT   = tab.log_ANT_RT_Alerting(tab_sel);


%% Inputs

pathway_contrasts = '/tedana009a1_vt/rsfc/';
contrast_names   = { 'ALFF_clean.nii'  'fALFF_clean.nii' } ;
groups.name = {'Control' 'Tinnitus'};
covars.name = {'Age','Genre','log_ANT_RT_Alerting'};

% getting scans & covariates organised per group
scans = cell(1,length(groups.name));
for igroup = 1:length(groups.name)
    j = 0 ;
    for iSubj = 1:length(e)
        
        ifile = e(iSubj).name;
        id = str2double(ifile(25:end)) ;
        subj_group = d.Groupe( find (d.Num_IRM == id) )  ;
        
        if subj_group == igroup
            j = j + 1 ;
            for icov = 1:length(covars.name)
                scans{igroup}.cov{j,1} = tab.Age(tab.code_IRM == id);
                scans{igroup}.cov{j,2} = tab.Genre(tab.code_IRM == id);
                scans{igroup}.cov{j,3} = tab.log_ANT_RT_Alerting(tab.code_IRM == id);
            end
            for icontr = 1:length(contrast_names) % missing RS folder path
                scans{igroup}.contrast{j, icontr } = fullfile( e(iSubj).getSerie('run_RS').path, pathway_contrasts, contrast_names{icontr}) ;
            end
            
        end
    end
end
groups.val = {scans{1}.contrast(:,1), scans{2}.contrast(:,1)}; % ALFF
covars.val = {vertcat(scans{1}.cov{:,1},scans{2}.cov{:,1});vertcat(scans{1}.cov{:,2},scans{2}.cov{:,2});vertcat(scans{1}.cov{:,3},scans{2}.cov{:,3})};

%% Model specify
clear par
par.sge = 0;
par.run = 1;
par.covars = 1;

varcov_2nd_level_2sample_model_spec(groups.val,outdirs(1),covars,par)

%% Model estimate
fspm = addsuffixtofilenames(outdirs(1),'/SPM.mat');
clear par
    par.run = 1;
    par.sge = 0;
    par.sge_queu = 'normal,bigmem';
    par.jobname  = 'spm_multireg_vbm_est';
    job_first_level_estimate(fspm,par)

%% Contrast definition
% F-stat
Main_effect_ANT = [0 0 0 1 0; 0 0 0 0 1];
Main_effect_ANT_ALFF = [0 0 0 1 -1];

contrast_F.names = {
    'Main_ANT_effect'
    'Main_ANT_effect_on_Group_fALFF'
    }';
contrast_F.values = {
    Main_effect_ANT
    Main_effect_ANT_ALFF
    }';
contrast_F.types = cat(1,repmat({'F'},[1,length(contrast_F.names)]));

contrast.names  = [contrast_F.names];
contrast.values = [contrast_F.values];
contrast.types  = [contrast_F.types];

%% Contrast : write
    clear par

    par.sge = 0;
    par.run = 1;
    par.display = 0;
    par.jobname = 'spm_write_vbm_con';

    par.sessrep = 'none';

    par.delete_previous = 1;
    par.report          = 0;
     
    job_first_level_contrast(fspm,contrast,par);



