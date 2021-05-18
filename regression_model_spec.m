function jobs = regression_model_spec(outdir, cons_a, cons_c, covars, target_regressor, par)
%%

  skip = [];
if ~exist ('par','var')
   par = '';
end
    
%% defpar
defpar.jobname                            = 'spm_spec_regression_auto';
defpar.walltime                           = '04:00:00';
defpar.run                                = 0;
defpar.nb_cond                            = 1;
defpar.nb_cons                            = 1;

par = complet_struct(par, defpar);
ijob = 1;
for i = 1 : par.nb_cond % a loop for each target regressor
    for iout = 1 : par.nb_cons % a loop for each contrast // con
      %jobs{ijob}.spm.stats.factorial_design.dir = {'/home/anna.skrzatek/data/nifti_test/resliced_ACT_gait/Rapid/rAPA_AP_duration/IL_V1'};
      jobs{ijob}.spm.stats.factorial_design.dir = outdir{i}(iout);

      jobs{ijob}.spm.stats.factorial_design.des.t2.scans1 = spm_select('expand',cons_a{iout}); 
      jobs{ijob}.spm.stats.factorial_design.des.t2.scans2 = spm_select('expand',cons_c{iout}); 
      jobs{ijob}.spm.stats.factorial_design.des.t2.dept = 0;
      jobs{ijob}.spm.stats.factorial_design.des.t2.variance = 1;
      jobs{ijob}.spm.stats.factorial_design.des.t2.gmsca = 0;
      jobs{ijob}.spm.stats.factorial_design.des.t2.ancova = 0;

      %% COVARIANTS - AGE(1) & GENDER(2)
      jobs{ijob}.spm.stats.factorial_design.cov(1).c = covars{1};
      %%
      jobs{ijob}.spm.stats.factorial_design.cov(1).cname = 'Age';
      jobs{ijob}.spm.stats.factorial_design.cov(1).iCFI = 1;
      jobs{ijob}.spm.stats.factorial_design.cov(1).iCC = 1;

      %%
      jobs{ijob}.spm.stats.factorial_design.cov(2).c = covars{2};
      %%
      jobs{ijob}.spm.stats.factorial_design.cov(2).cname = 'Gender';
      jobs{ijob}.spm.stats.factorial_design.cov(2).iCFI = 1;
      jobs{ijob}.spm.stats.factorial_design.cov(2).iCC = 1;

      %% TARGET REGRESSOR
      jobs{ijob}.spm.stats.factorial_design.cov(3).c = target_regressor.value{i};

      %jobs{ijob}.spm.stats.factorial_design.cov(3).cname = 'Delta_APA_AntPost';
      jobs{ijob}.spm.stats.factorial_design.cov(3).cname = target_regressor.name{i};

      jobs{ijob}.spm.stats.factorial_design.cov(3).iCFI = 2;
      jobs{ijob}.spm.stats.factorial_design.cov(3).iCC = 1;
      jobs{ijob}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
      jobs{ijob}.spm.stats.factorial_design.masking.tm.tm_none = 1;
      jobs{ijob}.spm.stats.factorial_design.masking.im = 1;
      jobs{ijob}.spm.stats.factorial_design.masking.em = {''};
      jobs{ijob}.spm.stats.factorial_design.globalc.g_omit = 1;
      jobs{ijob}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
      jobs{ijob}.spm.stats.factorial_design.globalm.glonorm = 1;

      spm('defaults','FMRI')
      ijob = ijob +1;
    %    job = do_cmd_sge(jobs{ijob},par);
    end
end

[ jobs ] = job_ending_rountines( jobs, skip, par );

end