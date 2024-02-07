%% paramètres / variables
clc
clear

main_dir = '/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/DATA/Non_chirurgicaux';
project_dir = '/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG';
cd (project_dir)
outdirs = {'/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Results/RS_2sample_ttest/ALFF','/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Results/RS_2sample_ttest/fALFF'} ;

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
            for icontr = 1:length(contrast_names)
                scans{igroup}.contrast{j, icontr } = fullfile( e(iSubj).path, pathway_contrasts, contrast_names{icontr}) ;
            end
            
        end
    end
end
groups.val = scans{1}.contrast(:,1) .toJob, scans{2}.contrast(:,1)}; % ALFF
covars.val = {vertcat(scans{1}.cov{:,1},scans{2}.cov{:,1});vertcat(scans{1}.cov{:,2},scans{2}.cov{:,2});vertcat(scans{1}.cov{:,3},scans{2}.cov{:,3})};

%% Model specify
clear par
par.sge = 0;
par.run = 1;
par.covars = 1;

varcov_2nd_level_2sample_model_spec(groups.val,outdirs(1),covars,par)