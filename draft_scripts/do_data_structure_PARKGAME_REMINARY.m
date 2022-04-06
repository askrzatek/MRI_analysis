function do_data_structure_PARKGAME_REMINARY(suj_source,root_desti)
%% INPUT
%% source subjects (nifti) folder
%% Destination folder, where all the subjects folder will be re-created with the symbolic links

% // either we load e.mat and we extract paths of files through examen architecture // may be easier to construct
% // either we create the architecture by getting each directory path by regex and at final step we create symbolic links to files from directories of which we extracted the path at step 1
% // step 1: extracting the paths for group separation --> PARK and REMINARY lists of paths
PARKGAME = get_subdir_regex(suj_source, 'PARKGAME');
REMINARY = get_subdir_regex(suj_source, 'REMINARY_\w{2}_');
% // step 2: extracting the paths for intra-group separation --> SESSION factor (V1/V2) lists of paths for each directory
P_V1 = get_subdir_regex(suj_source, 'PARKGAME.*1$');
P_V2 = get_subdir_regex(suj_source, 'PARKGAME.*2$');

R_V1 = get_subdir_regex('/home/anna.skrzatek/data/nifti/', 'REMINARY_\w{2}_.*1$');
R_V2 = get_subdir_regex('/home/anna.skrzatek/data/nifti/', 'REMINARY_\w{2}_.*2$');

%% // step 3: extracting the paths for models in each final folder of interest --> model_meica
model_PV1 = get_subdir_regex(P_V1, 'model_meica');
model_PV2 = get_subdir_regex(P_V2, 'model_meica');
cons_PV1 = get_subdir_regex_files(model_PV1,'^con');
cons_PV2 = get_subdir_regex_files(model_PV2,'^con');

model_RV1 = get_subdir_regex(R_V1, 'model_meica');
model_RV2 = get_subdir_regex(R_V2, 'model_meica');
cons_RV1 = get_subdir_regex_files(model_RV1,'^con');
cons_RV2 = get_subdir_regex_files(model_RV2,'^con');

%UNI = get_subdir_regex(suj_source, 'UNI_Images$');
%UNI = get_subdir_regex(suj_source, '^S09');
%UNI_im = get_subdir_regex_files(UNI,'.nii.gz');

%gre_field_map = get_subdir_regex(suj_source,'gre_field_mapping$');
%gre_field_map_phase = get_subdir_regex(suj_source,'gre_field_mapping_phase$');

%% On crée l'architecture
n = length(suj_source);


for k = 1:n
%     % Nom du sujet
    [path,~,ext] = fileparts(suj_source{k});
    [path,name,ext] = fileparts(path);
    suj = r_mkdir(root_desti, name);
    cmd='';
    
    % PARKGAME
%    cmd =  sprintf('%s\nln -f -s %s %s',cmd, PARKGAME{k}, suj{1});

    
    % REMINARY
%    cmd =  sprintf('%s\nln -f -s %s %s',cmd, REMINARY{k}, suj{1});

    
    % P_V1
    cmd =  sprintf('%s\nln -f -s %s %s',cmd, P_V1{k}, suj{1});
    
  
    % P_V2
    cmd =  sprintf('%s\nln -f -s %s %s',cmd, P_V2{k}, suj{1});

    
    % R_V1
    cmd =  sprintf('%s\nln -f -s %s %s',cmd, R_V1{k}, suj{1});

    
    % R_V2
    cmd =  sprintf('%s\nln -f -s %s %s',cmd, R_V2{k}, suj{1});

     
%    %INv1
%     cmd =  sprintf('%s\nln -f -s %s %s/%s/T1w_MP2RAGE_INV1.nii.gz',cmd, INV1_im{k},root_desti{1},name);
%     
%     % INv2
%     cmd =  sprintf('%s\nln -f -s %s %s/%s/T1w_MP2RAGE_INV2.nii.gz',cmd, INV2_im{k},root_desti{1},name);
%      
%     %UNI
%     cmd =  sprintf('%s\nln -f -s %s %s/%s/T1w_MP2RAGE_UNI.nii.gz',cmd, UNI_im{k},root_desti{1},name);
%     
%     if gre_field == 1
%     cmd =  sprintf('%s\ncp %s %s/%s',cmd, gre_field_map{k},root_desti{1},name);
%     cmd =  sprintf('%s\ncp %s %s/%s',cmd, gre_field_map_phase{k},root_desti{1},name);
%     end
%     
    unix(cmd);
end
end

