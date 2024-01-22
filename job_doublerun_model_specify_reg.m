function jobs = job_doublerun_model_specify_reg(dirFonc, dirOutput, onsets, emg_reg, par)

skip = [];
if ~exist ('par','var')
   par = '';
end


%% defpar
defpar.redo                               = 0;
defpar.jobname                            = 'spm_glm_spec_dble_auto';
defpar.walltime                           = '04:00:00';
defpar.run                                = 0;

par = complet_struct(par, defpar);

%%
for i = 1 : length(dirFonc)
    fspm = addsuffixtofilenames(dirOutput{i},'SPM.mat');
    if or(par.redo, ~exist(fspm{1},'file'))
        
        jobs{i}.spm.stats.fmri_spec.dir = dirOutput{i};
        jobs{i}.spm.stats.fmri_spec.timing.units = 'secs';
        jobs{i}.spm.stats.fmri_spec.timing.RT = 1.6;
        jobs{i}.spm.stats.fmri_spec.timing.fmri_t = 16;
        jobs{i}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
        %%
        subjectRun1 = get_subdir_regex_files(dirFonc{i,1},par.file_reg);
        unzip_volume(subjectRun1);
        subjectRun1 = get_subdir_regex_files(dirFonc{i,1}, par.file_reg, struct('verbose',0));

        clear allVolumes
        allVolumes = spm_select('expand',cellstr(subjectRun1{1}));
        jobs{i}.spm.stats.fmri_spec.sess(1).scans = allVolumes; %


        %%
        %jobs{i}.spm.stats.fmri_spec.sess(1).cond(1) = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
       if ~exist(emg_reg{i,1,1},'file') && ~exist(emg_reg{i,1,2},'file') && ~exist(emg_reg{i,1,3},'file') && ~exist(emg_reg{i,1,4},'file')
          % jobs{i}.spm.stats.fmri_spec.sess(1).multi = {'/home/anna.skrzatek/data/behav/2020_03_13_PARKGAMEII_044_CK_13_03_2020_V1_c/20200313T141827_PARKGAMEII_P044_CK_S01_13_03_2020_MRI_run01_SPM.mat'};
          jobs{i}.spm.stats.fmri_spec.sess(1).multi = onsets(i,1);
          jobs{i}.spm.stats.fmri_spec.sess(1).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
          jobs{i}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {}); %% maybe add the age & sex regressors ?

       else

          stim = load(onsets{i,1});
          
          %stim.names == "Rest"
          jobs{i}.spm.stats.fmri_spec.sess(1).cond(1).name = stim.names{1};
          jobs{i}.spm.stats.fmri_spec.sess(1).cond(1).onset = stim.onsets{1};
          jobs{i}.spm.stats.fmri_spec.sess(1).cond(1).duration = stim.durations{1};
          jobs{1}.spm.stats.fmri_spec.sess(1).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
          jobs{1}.spm.stats.fmri_spec.sess(1).cond(1).tmod = 0;
          jobs{1}.spm.stats.fmri_spec.sess(1).cond(1).orth = 1;
          
          %stim.names == "Imaginary_Left"
          jobs{i}.spm.stats.fmri_spec.sess(1).cond(2).name = stim.names{4};
          jobs{i}.spm.stats.fmri_spec.sess(1).cond(2).onset = stim.onsets{4};
          jobs{i}.spm.stats.fmri_spec.sess(1).cond(2).duration = stim.durations{4};
          jobs{1}.spm.stats.fmri_spec.sess(1).cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
          jobs{1}.spm.stats.fmri_spec.sess(1).cond(2).tmod = 0;
          jobs{1}.spm.stats.fmri_spec.sess(1).cond(2).orth = 1;

          %stim.names == "Imaginary_Right"
          jobs{i}.spm.stats.fmri_spec.sess(1).cond(3).name = stim.names{5};
          jobs{i}.spm.stats.fmri_spec.sess(1).cond(3).onset = stim.onsets{5};
          jobs{i}.spm.stats.fmri_spec.sess(1).cond(3).duration = stim.durations{5};
          jobs{1}.spm.stats.fmri_spec.sess(1).cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
          jobs{1}.spm.stats.fmri_spec.sess(1).cond(3).tmod = 0;
          jobs{1}.spm.stats.fmri_spec.sess(1).cond(3).orth = 1;


          %jobs{i}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {}); %% maybe add the age & sex regressors ?

          ext_R = load(emg_reg{i,1,1}); %% not a good indexation method // will depend on the structure I create in firstlevel_tedana_double_run_reg
          jobs{i}.spm.stats.fmri_spec.sess(1).regress(1).name = 'ext_R';
          jobs{i}.spm.stats.fmri_spec.sess(1).regress(1).val = ext_R.R;

          ext_L = load(emg_reg{i,1,2}); %% not a good indexation method
          jobs{i}.spm.stats.fmri_spec.sess(1).regress(2).name = 'ext_L';
          jobs{i}.spm.stats.fmri_spec.sess(1).regress(2).val = ext_L.R;

          fle_R = load(emg_reg{i,1,3});
          jobs{i}.spm.stats.fmri_spec.sess(1).regress(3).name = 'fle_R';
          jobs{i}.spm.stats.fmri_spec.sess(1).regress(3).val = fle_R.R;

          fle_L = load(emg_reg{i,1,4});
          jobs{i}.spm.stats.fmri_spec.sess(1).regress(4).name = 'fle_L';
          jobs{i}.spm.stats.fmri_spec.sess(1).regress(4).val = fle_L.R;

        end


        % jobs{i}.spm.stats.fmri_spec.sess(1).multi_reg = {'/home/anna.skrzatek/data/nifti_test/2020_03_13_PARKGAMEII_044_CK_13_03_2020_V1_c/S12_ACTIVATION/wts_multiple_regressors.txt'};
        if par.rp
            jobs{i}.spm.stats.fmri_spec.sess(1).multi_reg = get_subdir_regex_files(dirFonc{i,1},par.rp_regex);
        end
        jobs{i}.spm.stats.fmri_spec.sess(1).hpf = 128;


    %%
        subjectRun2 = get_subdir_regex_files(dirFonc{i,2},par.file_reg);
        unzip_volume(subjectRun2);
        subjectRun2 = get_subdir_regex_files(dirFonc{i,2}, par.file_reg, struct('verbose',0));
        clear allVolumes2
        allVolumes2 = spm_select('expand',cellstr(subjectRun2{1}));

        jobs{i}.spm.stats.fmri_spec.sess(2).scans = allVolumes2; %
        %%%%
        % % allVolumes = spm_select('expand',cellstr(subjectRuns{1}));
        % % jobs{1}.spm.stats.fmri_spec.sess(2).scans = allVolumes;
        %%
        %%
        %jobs{i}.spm.stats.fmri_spec.sess(1).cond(1) = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});

        if ~exist(emg_reg{i,2,1},'file')&& ~exist(emg_reg{i,2,2},'file') && ~exist(emg_reg{i,2,3},'file') && ~exist(emg_reg{i,2,4},'file')
          % jobs{i}.spm.stats.fmri_spec.sess(1).multi = {'/home/anna.skrzatek/data/behav/2020_03_13_PARKGAMEII_044_CK_13_03_2020_V1_c/20200313T141827_PARKGAMEII_P044_CK_S01_13_03_2020_MRI_run01_SPM.mat'};
          jobs{i}.spm.stats.fmri_spec.sess(2).multi = onsets(i,2);
          jobs{i}.spm.stats.fmri_spec.sess(2).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
          jobs{i}.spm.stats.fmri_spec.sess(2).regress = struct('name', {}, 'val', {}); %% maybe add the age & sex regressors ?

       else

          stim = load(onsets{i,2});

          %stim.names == "Rest"
          jobs{i}.spm.stats.fmri_spec.sess(2).cond(1).name = stim.names{1};
          jobs{i}.spm.stats.fmri_spec.sess(2).cond(1).onset = stim.onsets{1};
          jobs{i}.spm.stats.fmri_spec.sess(2).cond(1).duration = stim.durations{1};
          jobs{1}.spm.stats.fmri_spec.sess(2).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
          jobs{1}.spm.stats.fmri_spec.sess(2).cond(1).tmod = 0;
          jobs{1}.spm.stats.fmri_spec.sess(2).cond(1).orth = 1;
          
          %stim.names == "Imaginary_Left"
          jobs{i}.spm.stats.fmri_spec.sess(2).cond(2).name = stim.names{4};
          jobs{i}.spm.stats.fmri_spec.sess(2).cond(2).onset = stim.onsets{4};
          jobs{i}.spm.stats.fmri_spec.sess(2).cond(2).duration = stim.durations{4};
          jobs{1}.spm.stats.fmri_spec.sess(2).cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
          jobs{1}.spm.stats.fmri_spec.sess(2).cond(2).tmod = 0;
          jobs{1}.spm.stats.fmri_spec.sess(2).cond(2).orth = 1;

          %stim.names == "Imaginary_Right"
          jobs{i}.spm.stats.fmri_spec.sess(2).cond(3).name = stim.names{5};
          jobs{i}.spm.stats.fmri_spec.sess(2).cond(3).onset = stim.onsets{5};
          jobs{i}.spm.stats.fmri_spec.sess(2).cond(3).duration = stim.durations{5};
          jobs{1}.spm.stats.fmri_spec.sess(2).cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
          jobs{1}.spm.stats.fmri_spec.sess(2).cond(3).tmod = 0;
          jobs{1}.spm.stats.fmri_spec.sess(2).cond(3).orth = 1;

          %jobs{i}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {}); %% maybe add the age & sex regressors ?

          ext_R = load(emg_reg{i,2,1}); %% not a good indexation method // will depend on the structure I create in firstlevel_tedana_double_run_reg
          jobs{i}.spm.stats.fmri_spec.sess(2).regress(1).name = 'ext_R';
          jobs{i}.spm.stats.fmri_spec.sess(2).regress(1).val = ext_R.R;

          ext_L = load(emg_reg{i,2,2}); %% not a good indexation method
          jobs{i}.spm.stats.fmri_spec.sess(2).regress(2).name = 'ext_L';
          jobs{i}.spm.stats.fmri_spec.sess(2).regress(2).val = ext_L.R;

          fle_R = load(emg_reg{i,2,3});
          jobs{i}.spm.stats.fmri_spec.sess(2).regress(3).name = 'fle_R';
          jobs{i}.spm.stats.fmri_spec.sess(2).regress(3).val = fle_R.R;

          fle_L = load(emg_reg{i,2,4});
          jobs{i}.spm.stats.fmri_spec.sess(2).regress(4).name = 'fle_L';
          jobs{i}.spm.stats.fmri_spec.sess(2).regress(4).val = fle_L.R;

        end


        % jobs{i}.spm.stats.fmri_spec.sess(1).multi_reg = {'/home/anna.skrzatek/data/nifti_test/2020_03_13_PARKGAMEII_044_CK_13_03_2020_V1_c/S12_ACTIVATION/wts_multiple_regressors.txt'};
        if par.rp
            jobs{i}.spm.stats.fmri_spec.sess(2).multi_reg = get_subdir_regex_files(dirFonc{i,2},par.rp_regex);
        end
        jobs{i}.spm.stats.fmri_spec.sess(2).hpf = 128;

    %%

        jobs{i}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        jobs{i}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        jobs{i}.spm.stats.fmri_spec.volt = 1;
        jobs{i}.spm.stats.fmri_spec.global = 'None';
        jobs{i}.spm.stats.fmri_spec.mthresh = par.mask_thr;

        if par.mask
            if par.mask_path{1} == ''
                jobs{i}.spm.stats.fmri_spec.mask = unzip_volume(get_subdir_regex_files(dirFonc{i,2},'wbet.*mask'));
            else
                jobs{i}.spm.stats.fmri_spec.mask = unzip_volume(par.mask_path);
            end
        else
            jobs{i}.spm.stats.fmri_spec.mask = {''};
        end

        jobs{i}.spm.stats.fmri_spec.cvi = 'AR(1)';

        spm('defaults','FMRI')
    else
        skip = [skip i];
        fprintf('[%s]: skiping subj because %s exists \n',mfilename,fspm{1});
    end
%    job = do_cmd_sge(jobs{i},par);
end

[ jobs ] = job_ending_rountines( jobs, skip, par );

    
end