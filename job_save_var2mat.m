function jobs = test_batch_save()

jobs{1}.cfg_basicio.var_ops.cfg_save_vars.name = 'ROI';
jobs{1}.cfg_basicio.var_ops.cfg_save_vars.outdir = {'/network/iss/cenir/analyse/irm/users/anna.skrzatek/nifti/secondlevel_ACTIVATION_PARK_S1/ANOVA2x2_LxT/rois'};
jobs{1}.cfg_basicio.var_ops.cfg_save_vars.vars.vname = 'TabDat';
jobs{1}.cfg_basicio.var_ops.cfg_save_vars.vars.vcont = 'TabDat';
jobs{1}.cfg_basicio.var_ops.cfg_save_vars.saveasstruct = true;

spm_jobman('run', jobs);
end