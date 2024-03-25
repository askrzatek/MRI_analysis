function matlabbatch = varcov_2nd_level_2sample_model_spec(groups,outdirs,covars,par)
% TO BE TESTED %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Script for 2 sample t-test model specification : second level analysis
% Variable number of covariates & possible use of the param specifying if using 
% covariates & multiple covariates structure (whether we use external file or
% not - empty by default + variable number of covariates possible
%
% + paired version of the test (coming soon)
%
% A.Skrzatek Feb 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% groups : 2-cell structure where each cell contains strcells for T2 scans'
% paths per group
% outdirs : cell structure with string path to the output directory
% covars : structure with two subfields
% covars.name : strcell structure where each cell is the covariate name
% covars.val : cell structure (n,nsuj) where n is the number of covariates
% and nsuj the number of subjects data ordered per group (same order as in
% groups variable
% par : varying parameters for the model
% par.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
  skip = [];
if ~exist ('par','var')
   par = '';
end
    
%% defpar
defpar.jobname                            = 'spm_spec_secondlevel_RS_2s_auto';
defpar.walltime                           = '04:00:00';
defpar.run                                = 0;
defpar.paired                             = 0;
defpar.ancova                             = 0;
defpar.covars                             = 0;
if defpar.covars == 1
    defpar.intercov                           = cat(1,repmat(1,[1,length(covars.name)]));
end
defpar.multicov                           = struct('files', {}, 'iCFI', {}, 'iCC', {});

par = complet_struct(par, defpar);

%%
ijob = 1;
    
    jobs{ijob}.spm.stats.factorial_design.dir = outdirs(1);
    jobs{ijob}.spm.stats.factorial_design.des.t2.scans1 = spm_select('expand',groups{1});
    jobs{ijob}.spm.stats.factorial_design.des.t2.scans2 = spm_select('expand',groups{2});
    jobs{ijob}.spm.stats.factorial_design.des.t2.dept = 0;
    jobs{ijob}.spm.stats.factorial_design.des.t2.variance = 1;
    jobs{ijob}.spm.stats.factorial_design.des.t2.gmsca = 0;
    jobs{ijob}.spm.stats.factorial_design.des.t2.ancova = par.ancova;
%% Varying covariates loop : dependent on number of given covariates
    if par.covars == 0
        jobs{ijob}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    else
        for icov = 1:length(covars.name)
            jobs{ijob}.spm.stats.factorial_design.cov(icov).c = vertcat(covars.val{icov,:});
        %%
            jobs{ijob}.spm.stats.factorial_design.cov(icov).cname = covars.name{icov};
            jobs{ijob}.spm.stats.factorial_design.cov(icov).iCFI = par.intercov(icov);
            % par.intercov(icov) == 
            %                       1: no interaction 
            %                       2: interaction with factor 1
            %                       3: interaction with factor 2
            %                       4: interaction with factor 3
            jobs{ijob}.spm.stats.factorial_design.cov(icov).iCC = 1;
            %%
        end
    end
    jobs{ijob}.spm.stats.factorial_design.multi_cov = par.multicov;
    jobs{ijob}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    jobs{ijob}.spm.stats.factorial_design.masking.im = 1;
    jobs{ijob}.spm.stats.factorial_design.masking.em = {''};
    jobs{ijob}.spm.stats.factorial_design.globalc.g_omit = 1;
    jobs{ijob}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    jobs{ijob}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    spm('defaults','FMRI')
    

[ jobs ] = job_ending_rountines( jobs, skip, par );

end