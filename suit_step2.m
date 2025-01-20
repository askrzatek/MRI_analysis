clearvars;
clc;

suit_cerebellar_lobules = '/network/iss/cenir/analyse/irm/users/cecile.gallea/ASYA/asyasuit/spm12/toolbox/suit/atlasesSUIT/Lobules-SUIT.nii';
results_path = '/network/iss/cenir/analyse/irm/users/anna.skrzatek/nifti_test/PARKGAMEIIsuit_results';
if ~exist(results_path, 'dir'), mkdir(results_path), end

cerebellum_stats = fullfile(results_path, 'suit_cerebellar_lobules.csv');

if ~exist(cerebellum_stats, 'file')
    cmd = '';
    cmd = sprintf('%s\necho "" >> %s\n', cmd, cerebellum_stats); 
    cmd = sprintf('%s\necho -n subjectID, Left_I_IV, Right_I_IV, Left_V, Right_V, Left_VI, Vermis_VI, Right_VI, Left_CrusI, Vermis_CrusI, Right_CrusI, Left_CrusII, >> %s\n', cmd, cerebellum_stats);
    cmd = sprintf('%s\necho -n Vermis_CrusII, Right_CrusII, Left_VIIb, Vermis_VIIb, Right_VIIb, Left_VIIIa, Vermis_VIIIa, Right_VIIIa, Left_VIIIb, Vermis_VIIIb, Right_VIIIb, >> %s\n', cmd, cerebellum_stats);
    cmd = sprintf('%s\necho -n Left_IX, Vermis_IX, Right_IX, Left_X, Vermis_X, Right_X, Left_Dentate, Right_Dentate, Left_Interposed, Right_Interposed, Left_Fastigial, Right_Fastigial >> %s\n', cmd, cerebellum_stats);
    unix(cmd);
end


suj_path = '/network/iss/cenir/analyse/irm/users/anna.skrzatek/nifti_test/';

suj = get_subdir_regex(suj_path);

[pp, name] = get_parent_path(suj);
par.verbose = 0; 


for i = 1 : length(suj)
    temp = regexp(name{i}, '_', 'split');
    sujID = [temp{end-6}, '_', temp{end-5}, '_', temp{end-1}, '_', temp{end}];
    
    %% suit isolation into gray, white matter and pcereb (7 mins)
    suitdir = fullfile(suj{i}, 'suit');
    t1_file = char(get_subdir_regex_files(suitdir, ['^', sujID, '_t1corr.nii'], par));
    t1_file = remove_gz(t1_file); 
    
    if isempty(t1_file)
        fprintf('%s has no t1. skipping...\n', sujID);
        continue
    end
    
    if exist(fullfile(suitdir, ['c_', sujID, '_t1corr_pcereb.nii']), 'file')
        fprintf('--%s: segmentation was previously done\n', sujID)
    else

        t2_file = char(get_subdir_regex_files(suitdir, 't2regWarped', par));
        if ~isempty(t2_file)
            fprintf('--%s : using t1 and t2 files in suit_isolation_seg\n', sujID);
            t2_file = remove_gz(t2_file);
            suit_isolate_seg({t1_file, t2_file}, 'maskp', 0.1);
        else
            fprintf('--%s : using only t1 file in suit_isolation_seg\n', sujID);
            suit_isolate_seg({t1_file}, 'maskp', 0.1);
        end

        %This pause is included to allow for verification of pcereb mask
        fprintf('\n######please check and edit the pcereb for %s\n', sujID);
        fprintf('---hit enter when done\n');
        pause;
    end
    
    %% suite normalization with denate roi (3 mins)
    dentate_roi = char(get_subdir_regex_files(suitdir, 'dentate_roi', par));
    qsm_reg = char(get_subdir_regex_files(suitdir, 'qsmreg.nii', par)); 
    
    if isempty(dentate_roi) || isempty(qsm_reg)
        fprintf('%s has no dentate_roi. normalizing with no dentate...\n', sujID);
        
        % normalization without dentate mask
        job = '';
        job.subjND.gray = {char(get_subdir_regex_files(suitdir, 'seg1', par))};
        job.subjND.white = {char(get_subdir_regex_files(suitdir, 'seg2', par))};
        job.subjND.isolation = {char(get_subdir_regex_files(suitdir, 'pcereb', par))};

        suit_normalize_dartel(job);  
                 
        % reslice to subject space
        fprintf('--%s : reslicing to native space\n', sujID);
        job = '';
        job.Affine = {char(get_subdir_regex_files(suitdir, '^Affine', par))};
        job.flowfield = {char(get_subdir_regex_files(suitdir, '^u_a', par))};
        job.resample = {suit_cerebellar_lobules};
        job.ref = {t1_file};
        suit_reslice_dartel_inv(job);
        
        
        % Get individual lobular volume
        native_lobules = char(get_subdir_regex_files(suitdir, 'iw_Lobules', par));
        
        cmd = '';
        cmd = sprintf('%s\necho "" >> %s\n', cmd, cerebellum_stats);
        cmd = sprintf('%s\necho -n %s, >> %s\n', cmd, sujID, cerebellum_stats);
        
        l_thresh = 0; u_thresh = 2;
        for kk = 1 : 34
            vol_cmd = sprintf('fslstats %s -l %d -u %d -V \n', native_lobules, l_thresh, u_thresh);
            [~, lob_vol] = unix(vol_cmd);
            lob_vol = regexp(lob_vol, ' ', 'split');
            lob_vol = lob_vol{1};

            cmd = sprintf('%s\necho -n %s, >> %s \n', cmd, lob_vol, cerebellum_stats);            
            l_thresh = l_thresh + 1;
            u_thresh = u_thresh + 1;
            
        end
        
        unix(cmd);
        
        
    else
        fprintf('--%s : normalization with dentate roi\n', sujID);
        dentate_roi = remove_gz(dentate_roi);

        job = '';
        job.subjND.gray = {char(get_subdir_regex_files(suitdir, 'seg1', par))};
        job.subjND.white = {char(get_subdir_regex_files(suitdir, 'seg2', par))};
        job.subjND.isolation = {char(get_subdir_regex_files(suitdir, 'pcereb', par))};
        job.subjND.dentateROI = {dentate_roi};

        suit_normalize_dentate(job);    
        
        % reslice to subject space
        fprintf('--%s : reslicing to native space\n', sujID);
        job = '';
        job.Affine = {char(get_subdir_regex_files(suitdir, '^Affine', par))};
        job.flowfield = {char(get_subdir_regex_files(suitdir, '^u_a', par))};
        job.resample = {suit_cerebellar_lobules};
        job.ref = {t1_file};
        suit_reslice_dartel_inv(job);
        
        % Get whole cerebellar volume
        native_cerebellum = char(get_subdir_regex_files(suitdir, 'iw_Cerebellum', par));
        
        cmd = '';
        cmd = sprintf('%s\necho "" >> %s\n', cmd, cerebellum_stats);
        cmd = sprintf('%s\necho -n %s, >> %s\n', cmd, sujID, cerebellum_stats);
        cmd = sprintf('%s\necho -n " " >> %s\n', cmd, cerebellum_stats);
        
        vol_cmd = sprintf('fslstats %s -V \n', native_cerebellum);
        [~, cereb_vol] = unix(vol_cmd);
        cereb_vol = regexp(cereb_vol, ' ', 'split');
        cereb_vol = cereb_vol{1};
        cmd = sprintf('%s\necho -n %s, >> %s \n', cmd, cereb_vol, cerebellum_stats);
        cmd = sprintf('%s\necho -n " " >> %s\n', cmd, cerebellum_stats);
        
        
        % Get individual lobular volume
        native_lobules = char(get_subdir_regex_files(suitdir, 'iw_Lobules', par));
        
        l_thresh = 0; u_thresh = 2;
        for kk = 1 : 34
            vol_cmd = sprintf('fslstats %s -l %d -u %d -V \n', native_lobules, l_thresh, u_thresh);
            [~, lob_vol] = unix(vol_cmd);
            lob_vol = regexp(lob_vol, ' ', 'split');
            lob_vol = lob_vol{1};
            cmd = sprintf('%s\necho -n %s, >> %s \n', cmd, lob_vol, cerebellum_stats);
            cmd = sprintf('%s\necho -n " " >> %s\n', cmd, cerebellum_stats);
            
            l_thresh = l_thresh + 1;
            u_thresh = u_thresh + 1;
            
        end
        
        unix(cmd);

        %### reslice to template space (3 mins)
% %         job = '';
% %         job.subj.affineTr = {char(get_subdir_regex_files(suitdir, '^Affine', par))};
% %         job.subj.flowfield = {char(get_subdir_regex_files(suitdir, '^u_a', par))};
% %         job.subj.resample = {t1_file};
% %         job.subj.mask = {char(get_subdir_regex_files(suitdir, 'pcereb', par))};
% %         suit_reslice_dartel(job); %whole cerebellum
        
% %         %summary of t1 cerebellar regions
% %         t1_to_suit = char(get_subdir_regex_files(suitdir, ['^wd.*t1corr',sujID], par));
% %         t1_summary = fullfile(suitdir, [sujID, '_t1summary.txt']);
% %         suit_ROI_summarize(t1_to_suit, 'atlas', suit_cerebellar_lobules, 'stats', {'size'}, 'outfilename', t1_summary);
% %           
        % reslice dentate to template space
        qsm_reg = remove_gz(qsm_reg);
        
        job = '';
        job.subj.affineTr = {char(get_subdir_regex_files(suitdir, '^Affine', par))};
        job.subj.flowfield = {char(get_subdir_regex_files(suitdir, '^u_a', par))};
        job.subj.resample = {qsm_reg};
        job.subj.mask = {dentate_roi};
        suit_reslice_dartel(job);
        
        dentate_to_suit = char(get_subdir_regex_files(suitdir, '^wd.*qsmreg', par));
        
        check_dim = sprintf('fslval %s dim1\n', dentate_to_suit);
        [~, xdim] = unix(check_dim);
        xdim = str2double(xdim);
        
        %check for the mid sagittal slice to separate right and left
        %dentate
        mid_sagittal_slice = round(xdim/2,0);
        
        fprintf('--%s : extracting cerebellar volumes and mean qsm in dentate\n', sujID);
        
        left_dentate = fullfile(suitdir, [sujID, '_left_dentate.nii.gz']);
        right_dentate = fullfile(suitdir, [sujID, '_right_dentate.nii.gz']);
                
        cmd = '';
        cmd = sprintf('%s\nfslmaths %s -roi 0 %d 0 -1 0 -1 0 1 %s\n', cmd, dentate_to_suit, mid_sagittal_slice, left_dentate);
        cmd = sprintf('%s\nfslmaths %s -roi %d -1 0 -1 0 -1 0 1 %s\n', cmd, dentate_to_suit, mid_sagittal_slice, right_dentate);
        unix(cmd);
        
        % volume of dentate. lower threshold = 5
        cmd = '';
        cmd = sprintf('%s\necho "" >> %s\n', cmd, dentate_stats_output);
        cmd = sprintf('%s\necho -n %s, >> %s\n', cmd, sujID, dentate_stats_output);
        cmd = sprintf('%s\necho -n " " >> %s\n', cmd, dentate_stats_output);
        
        vol_cmd = sprintf('fslstats %s -l 5 -V \n', left_dentate);
        [~, lvolume] = unix(vol_cmd);
        lvolume = regexp(lvolume, ' ', 'split');
        lvolume = lvolume{1};
        cmd = sprintf('%s\necho -n %s, >> %s \n', cmd, lvolume, dentate_stats_output);
        cmd = sprintf('%s\necho -n " " >> %s\n', cmd, dentate_stats_output);
        
        vol_cmd = sprintf('fslstats %s -l 5 -V \n', right_dentate);
        [~, rvolume] = unix(vol_cmd);
        rvolume = regexp(rvolume, ' ', 'split');
        rvolume = rvolume{1};
        cmd = sprintf('%s\necho -n %s, >> %s \n', cmd, rvolume, dentate_stats_output);
        
        % mean qsm. lower threshold = 33
        cmd = sprintf('%s\necho -n " " >> %s\n', cmd, dentate_stats_output);
        cmd = sprintf('%s\nmean=$(fslstats %s -l 33 -M) \n', cmd, left_dentate);
        cmd = sprintf('%s\necho -n "$mean," >> %s \n', cmd, dentate_stats_output);
        cmd = sprintf('%s\necho -n " " >> %s\n', cmd, dentate_stats_output);
        cmd = sprintf('%s\nmean=$(fslstats %s -l 33 -M) \n', cmd, right_dentate);
        cmd = sprintf('%s\necho -n "$mean," >> %s \n', cmd, dentate_stats_output);
        
        % histograms
        cmd = sprintf('%s\necho "" >> %s\n', cmd, left_dentate_hist);
        cmd = sprintf('%s\necho -n %s, >> %s\n', cmd, sujID, left_dentate_hist);
        hist_cmd = sprintf('fslstats %s -l 33 -h 35 \n', left_dentate);
        [~, hist_output] = unix(hist_cmd);
        hist_output = regexprep(hist_output,'[\n\r\ ]+',',');
        cmd = sprintf('%s\necho -n %s >> %s \n', cmd, hist_output, left_dentate_hist);
        
        cmd = sprintf('%s\necho "" >> %s\n', cmd, right_dentate_hist);
        cmd = sprintf('%s\necho -n %s, >> %s\n', cmd, sujID, right_dentate_hist);
        hist_cmd = sprintf('fslstats %s -l 33 -h 35 \n', right_dentate);
        [~, hist_output] = unix(hist_cmd);
        hist_output = regexprep(hist_output,'[\n\r\ ]+',',');
        cmd = sprintf('%s\necho -n %s >> %s \n', cmd, hist_output, right_dentate_hist);

        unix(cmd);
        
    end
    
    if i == length(suj)
        fprintf('\n=====all done, enjoy!!========\n');
    end
    
    
end

