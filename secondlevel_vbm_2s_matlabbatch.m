function matlabbatch = secondlevel_vbm_2s_matlabbatch(groups,outdirs,covars,par)
%%

  skip = [];
if ~exist ('par','var')
   par = '';
end
    
%% defpar
defpar.jobname                            = 'spm_spec_secondlevel_RS_2s_auto';
defpar.walltime                           = '04:00:00';
defpar.run                                = 0;

par = complet_struct(par, defpar);

%%
ijob = 1;
    
    jobs{ijob}.spm.stats.factorial_design.dir = outdirs{1};
    jobs{ijob}.spm.stats.factorial_design.des.t2.scans1 = spm_select('expand',groups{1});
    jobs{ijob}.spm.stats.factorial_design.des.t2.scans2 = spm_select('expand',groups{2});
    jobs{ijob}.spm.stats.factorial_design.des.t2.dept = 0;
    jobs{ijob}.spm.stats.factorial_design.des.t2.variance = 1;
    jobs{ijob}.spm.stats.factorial_design.des.t2.gmsca = 0;
    jobs{ijob}.spm.stats.factorial_design.des.t2.ancova = 0;
    jobs{ijob}.spm.stats.factorial_design.cov(1).c = vertcat(covars{1,:});
    jobs{ijob}.spm.stats.factorial_design.cov(1).cname = 'Age';
    jobs{ijob}.spm.stats.factorial_design.cov(1).iCFI = 1;
    jobs{ijob}.spm.stats.factorial_design.cov(1).iCC = 1;
    jobs{ijob}.spm.stats.factorial_design.cov(2).c = vertcat(covars{2,:});
    jobs{ijob}.spm.stats.factorial_design.cov(2).cname = 'Gender';
    jobs{ijob}.spm.stats.factorial_design.cov(2).iCFI = 1;
    jobs{ijob}.spm.stats.factorial_design.cov(2).iCC = 1;
    jobs{ijob}.spm.stats.factorial_design.cov(3).c = vertcat(covars{1,:});
    jobs{ijob}.spm.stats.factorial_design.cov(3).cname = 'TIV';
    jobs{ijob}.spm.stats.factorial_design.cov(3).iCFI = 1;
    jobs{ijob}.spm.stats.factorial_design.cov(3).iCC = 1;

    %     jobs{ijob}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    jobs{ijob}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    jobs{ijob}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    jobs{ijob}.spm.stats.factorial_design.masking.im = 1;
    jobs{ijob}.spm.stats.factorial_design.masking.em = {''};
    jobs{ijob}.spm.stats.factorial_design.globalc.g_omit = 1;
    jobs{ijob}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    jobs{ijob}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    spm('defaults','FMRI')
    

[ jobs ] = job_ending_rountines( jobs, skip, par );

end