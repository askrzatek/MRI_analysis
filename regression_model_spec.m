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

par = complet_struct(par, defpar);

for i = 1 : par.nb_cond
  %jobs{i}.spm.stats.factorial_design.dir = {'/home/anna.skrzatek/data/nifti_test/resliced_ACT_gait/Rapid/rAPA_AP_duration/IL_V1'};
  jobs{i}.spm.stats.factorial_design.dir = outdir{i};

  jobs{i}.spm.stats.factorial_design.des.t2.scans1 = spm_select('expand',cons_a); % needs testing : whether a loop is necessary or we can get by
  jobs{i}.spm.stats.factorial_design.des.t2.scans2 = spm_select('expand',cons_c); % needs testing : whether a loop is necessary or we can get by
  jobs{i}.spm.stats.factorial_design.des.t2.dept = 0;
  jobs{i}.spm.stats.factorial_design.des.t2.variance = 1;
  jobs{i}.spm.stats.factorial_design.des.t2.gmsca = 0;
  jobs{i}.spm.stats.factorial_design.des.t2.ancova = 0;
  
  %% COVARIANTS - AGE(1) & GENDER(2)
  %jobs{i}.spm.stats.factorial_design.cov(1).c = [70
  %                                                      74
  %                                                      64
  %                                                      76
  %                                                      61
  %                                                      75
  %                                                      66
  %                                                      72
  %                                                      68
  %                                                      68
  %                                                      72
  %                                                      56
  %                                                      57];
  jobs{i}.spm.stats.factorial_design.cov(1).c = covars{1}
  %%
  jobs{i}.spm.stats.factorial_design.cov(1).cname = 'Age';
  jobs{i}.spm.stats.factorial_design.cov(1).iCFI = 1;
  jobs{i}.spm.stats.factorial_design.cov(1).iCC = 1;
  
  %%
  %jobs{i}.spm.stats.factorial_design.cov(2).c = [1
  %                                                      1
  %                                                      1
  %                                                      2
  %                                                      1
  %                                                      2
  %                                                      2
  %                                                      1
  %                                                      2
  %                                                      2
  %                                                      2
  %                                                      2
  %                                                      1];
  jobs{i}.spm.stats.factorial_design.cov(2).c = covars{2}
  %%
  jobs{i}.spm.stats.factorial_design.cov(2).cname = 'Gender';
  jobs{i}.spm.stats.factorial_design.cov(2).iCFI = 1;
  jobs{i}.spm.stats.factorial_design.cov(2).iCC = 1;
  
  %% TARGET REGRESSOR
  % jobs{i}.spm.stats.factorial_design.cov(3).c = [2.39874527879115
  %                                                      -7.76535427782036
  %                                                      5.7803799483553
  %                                                      18.0505768971419
  %                                                      5.36499247489797
  %                                                      14.7903939438382
  %                                                      3.54440383781971
  %                                                      6.35406411361517
  %                                                      12.0672220288013
  %                                                      12.3199655881826
  %                                                      -2.70590490795713
  %                                                      -15.8568741963904
  %                                                      11.3793264818906];
  jobs{i}.spm.stats.factorial_design.cov(3).c = target_regressor.value;
  
  %jobs{i}.spm.stats.factorial_design.cov(3).cname = 'Delta_APA_AntPost';
  jobs{i}.spm.stats.factorial_design.cov(3).cname = target_regressor.name;

  jobs{i}.spm.stats.factorial_design.cov(3).iCFI = 2;
  jobs{i}.spm.stats.factorial_design.cov(3).iCC = 1;
  jobs{i}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
  jobs{i}.spm.stats.factorial_design.masking.tm.tm_none = 1;
  jobs{i}.spm.stats.factorial_design.masking.im = 1;
  jobs{i}.spm.stats.factorial_design.masking.em = {''};
  jobs{i}.spm.stats.factorial_design.globalc.g_omit = 1;
  jobs{i}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
  jobs{i}.spm.stats.factorial_design.globalm.glonorm = 1;

  spm('defaults','FMRI')

%    job = do_cmd_sge(jobs{i},par);
end

[ jobs ] = job_ending_rountines( jobs, skip, par );

end