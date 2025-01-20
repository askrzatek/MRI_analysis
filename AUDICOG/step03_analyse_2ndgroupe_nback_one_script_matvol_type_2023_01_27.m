%% paramètres / variables
clc
clear

main_dir = '/network/iss/cenir/analyse/irm/studies/AUDICOG/DATA/Non_chirurgicaux';
cd '/network/iss/cenir/analyse/irm/studies/AUDICOG';

load('e.mat');

% fichier de correspondance numero IRM - comportement - groupe - age
d = readtable( [ './DATA/' , 'Correspondance_Numero_Comportement_IRM.csv' ])  ;

output_dir = '/network/iss/cenir/analyse/irm/studies/AUDICOG/Results/NBack' ;

%% Pathway

pathway_contrasts = '/stats/nback_01/';
contrast_names   = { 'con_0004.nii'  'con_0005.nii' } ;   % Correspond 0004 : 0back / 0005 2back
conditions.names = {'Groupe' 'nBack' } ;
conditions.levels  = [ 2  2] ; % Nombre de niveau pour chaque condition

% Prepare the list of scans
for igroup = 1:conditions.levels(1)
    j = 0 ;
    for iSubj = 1:length(e)
        
        ifile = e(iSubj).name;
        id = str2double(ifile(25:end)) ;
        subj_group = d.Groupe( find (d.Num_IRM == id) )  ;
        
        if subj_group == igroup
            j = j + 1 ;
            for icontr = 1:conditions.levels(2)
                scans(igroup).contrast{j, icontr } = fullfile( e(iSubj).path, pathway_contrasts, contrast_names{icontr}) ;
            end
            
        end
    end
end

%% Prepare the matlabbatch
k = 0;
for igroup = 1:conditions.levels(1)
    for icontr = 1:conditions.levels(2)
        k = k + 1 ;
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(k).levels = [igroup icontr];
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(k).scans =  scans(igroup).contrast(: ,icontr ) ;
        
    end
end

matlabbatch{1}.spm.stats.factorial_design.dir = {output_dir};

% Info sur le facteur 1 = Groupe
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = conditions.names{1};
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = conditions.levels(1);
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;

% Info sur le facteur 2 = Frequences auditives
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name =conditions.names{2};
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).levels = conditions.levels(2);
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).ancova = 0;

%
matlabbatch{1}.spm.stats.factorial_design.des.fd.contrasts = 1;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;


clear par
par.run = 1;
par.sge = 0;
par.display = 0;

skip = [] ;

job_ending_rountines(matlabbatch, skip, par);



%% Estimate

clear par
par.run = 1;
par.sge = 0;

spmmat = fullfile(output_dir,'SPM.mat');

matlabbatch{1}.spm.stats.fmri_est.spmmat = {spmmat} ;
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;


skip = [] ;

job_ending_rountines(matlabbatch,skip,par);


