function matlabbatch = varcov_2nd_level_ANOVA_specify(factors, outdirs, covars, par)
% TO BE TESTED %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Script for ANOVA model specification : second level analysis
% Variable number of factors & levels + possible use of the param specifying if
% using covariates & the multiple covariates structure (whether we use external
% file or not) + variable number of covariates possible
%
% A.Skrzatek Feb 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% outdirs : cell structure with string path to the output directory
% factors : structure with two subfields
% factors.name : strcell structure where each cell is the factor name
% factors.val : cell structure (nfactor x nlevel) where nfactor is the number of
% factors and nlevel the number of levels per factor each cell contains
% strcell paths to T2 scans per condition, indexes correspond to condition
% ID : factors.val{1,1} is the factor 1, level 1 scan list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check input arguments
    skip = [];
if ~exist('par','var')
    par = ''; % for defpar
end

%% defpar
defpar.jobname                                = 'spm_ANOVA_model_spec';
defpar.walltime                               = '04:00:00';
defpar.run                                    = 0;
defpar.nfactor                                = 2;
defpar.nlevels                                = [2 2];
defpar.ancova                                 = [0 0];
defpar.covars                                 = 0;
defpar.multicov                               = struct('files', {}, 'iCFI', {}, 'iCC', {});

par = complet_struct(par, defpar);

%%
ijob = 1;

    jobs{ijob}.spm.stats.factorial_design.dir = outdirs{ijob};
    %% Level definition and corresponding scans loading from 'factors.name' structure
    for ifactor = 1: length(factors.name)
        jobs{ijob}.spm.stats.factorial_design.des.fd.fact(ifactor).name = factors.name{ifactor};
        jobs{ijob}.spm.stats.factorial_design.des.fd.fact(ifactor).levels = par.nlevel(ifactor);
        jobs{ijob}.spm.stats.factorial_design.des.fd.fact(ifactor).dept = 0;
        jobs{ijob}.spm.stats.factorial_design.des.fd.fact(ifactor).variance = 1;
        jobs{ijob}.spm.stats.factorial_design.des.fd.fact(ifactor).gmsca = 0;
        jobs{ijob}.spm.stats.factorial_design.des.fd.fact(ifactor).ancova = par.ancova(ifactor);
    end
    %% Level definition and corresponding scans loading from 'factors.val' structure
    lcell = 1;
    for ilevel = 1: par.nfactor
        for ilevel2 = 1: par.nlevel(ilevel)
            jobs{ijob}.spm.stats.factorial_design.des.fd.icell(lcell).levels = [ilevel
                                                                                ilevel2];
            jobs{ijob}.spm.stats.factorial_design.des.fd.icell(lcell).scans = factors.val{ilevel,ilevel2};
            lcell = l+1;
        end
    end
    
    jobs{ijob}.spm.stats.factorial_design.des.fd.contrasts = 0;
    %% Varying covariates loop : dependent on number of given covariates
    if par.covars == 0
        jobs{ijob}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    else
        for icov = 1:length(covars.name)
            jobs{ijob}.spm.stats.factorial_design.cov(icov).c = vertcat(covars.val{icov,:});
        %%
            jobs{ijob}.spm.stats.factorial_design.cov(icov).cname = covars.name{icov};
            jobs{ijob}.spm.stats.factorial_design.cov(icov).iCFI = 1;
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
%spm_jobman('run',jobs)

end
