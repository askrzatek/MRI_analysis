function jobs = job_doublerun_model_specify_reg(dirFonc, dirOutput, onsets, emg_reg, par)

skip = [];
if ~exist ('par','var')
   par = '';
end


%% defpar
defpar.jobname                            = 'spm_glm_spec_dble_auto';
defpar.walltime                           = '04:00:00';
defpar.run                                = 0;

par = complet_struct(par, defpar);

%%
for i = 1 : length(dirFonc)
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
   if ~exist emg_reg
      % jobs{i}.spm.stats.fmri_spec.sess(1).multi = {'/home/anna.skrzatek/data/behav/2020_03_13_PARKGAMEII_044_CK_13_03_2020_V1_c/20200313T141827_PARKGAMEII_P044_CK_S01_13_03_2020_MRI_run01_SPM.mat'};
      jobs{i}.spm.stats.fmri_spec.sess(1).multi = onsets(i,1);
   else


      %% QUESTION : comment organiser les emg_reg comme variable - 2 fichiers par session par sujet
      %cond1 = load(emg_reg(i,1));

      cond1 = load(emg_reg{i,1}(1))
      jobs{i}.spm.stats.fmri_spec.sess(1).cond(1).onset = cond1.onset; % or maybe cond1(:,1) ?

      %jobs{i}.spm.stats.fmri_spec.sess(1).cond(2) = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
      %cond2 = load(emg_reg(i,2));
      cond2 = load(emg_reg{i,1}(2))
      jobs{i}.spm.stats.fmri_spec.sess(1).cond(2).onset = cond2.onset; % or maybe cond2(:,1) ?
    %% END OF QUESTION

    
      %% QUESTION : utiliser comme multi ou rajouter les conditions if positives pour creer les 2 conditions qui manquent ?
      
      stim = load(onsets(i,1));

      if stim.name == ""
         jobs{i}.spm.stats.fmri_spec.sess(1).cond(3).onset = stim.onset;
      end

      if stim.name == ""
         jobs{i}.spm.stats.fmri_spec.sess(1).cond(4).onset = stim.onset;
      end
    %% END OF QUESTION

    end

    jobs{i}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {}); %% maybe add the age & sex regressors ?

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

    if ~exist emg_reg
         % jobs{i}.spm.stats.fmri_spec.sess(2).multi = {'/home/anna.skrzatek/data/behav/2020_03_13_PARKGAMEII_044_CK_13_03_2020_V1_c/20200313T141827_PARKGAMEII_P044_CK_S01_13_03_2020_MRI_run01_SPM.mat'};
         jobs{i}.spm.stats.fmri_spec.sess(2).multi = onsets(i,2);
    else

       %% QUESTION : comment organiser les emg_reg comme variable - 2 fichiers par session par sujet
       cond1 = load(emg_reg{i,2}(1));
       jobs{i}.spm.stats.fmri_spec.sess(2).cond(1).onset = cond1.onset;

       %jobs{i}.spm.stats.fmri_spec.sess(1).cond(2) = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
       cond2 = load(emg_reg{i,2}(2));
       jobs{i}.spm.stats.fmri_spec.sess(2).cond(2).onset = cond2.onset;
       %% END OF QUESTION

       
       %% QUESTION : utiliser comme multi ou rajouter les conditions if positives pour creer les 2 conditions qui manquent ?
       % jobs{i}.spm.stats.fmri_spec.sess(1).multi = {'/home/anna.skrzatek/data/behav/2020_03_13_PARKGAMEII_044_CK_13_03_2020_V1_c/20200313T141827_PARKGAMEII_P044_CK_S01_13_03_2020_MRI_run01_SPM.mat'};
       stim = load(onsets(i,2));

       if stim.name == ""
         jobs{i}.spm.stats.fmri_spec.sess(2).cond(3).onset = stim.onset;
       else skip = [skip stim.name]
       end

       if stim.name == ""
         jobs{i}.spm.stats.fmri_spec.sess(2).cond(4).onset = stim.onset;
       else skip = [skip stim.name]
       end
       %% END OF QUESTION

    end

    jobs{i}.spm.stats.fmri_spec.sess(2).regress = struct('name', {}, 'val', {}); %% maybe add the age & sex regressors ?

    % jobs{i}.spm.stats.fmri_spec.sess(1).multi_reg = {'/home/anna.skrzatek/data/nifti_test/2020_03_13_PARKGAMEII_044_CK_13_03_2020_V1_c/S12_ACTIVATION/wts_multiple_regressors.txt'};
    if par.rp
        jobs{i}.spm.stats.fmri_spec.sess(2).multi_reg = get_subdir_regex_files(dirFonc{i,1},par.rp_regex);
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

%    job = do_cmd_sge(jobs{i},par);
end

[ jobs ] = job_ending_rountines( jobs, skip, par );

    
end