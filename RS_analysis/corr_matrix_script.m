
suj = getSubjects;

dmodel = gdir(suj,'^GLM_NativeSapce') 

%dmodel = gdir(suj,'GLM_RSnormalised')   
%glm = gdir(dmodel,'GLM2_s4'); 

glm = gdir(dmodel,'GLM2_s6'); 




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% voi = get_subdir_regex_files(glm,'VOI.*mat'); 
% 
% [diff,matched] = compare_cellnames2ref(voi);
% 
% fmat = get_subdir_regex_files(glm,matched); 
% 
% idr = {'Right_';'Right-'; '_rh'; '_RH_';'right-';'right_';'-rh-'; '_R_'}
% idl = {'Left_';'Left-'; '_lh'; '_LH_';'left-';'left_'; '-lh-'; '_L_'}
% 
% [notInIDr,notInfmatr, matchr] = compare_suj_cells(cellstr(fmat{2}),idr) % region without overlap
% [notInIDl,notInfmatl, matchl] = compare_suj_cells(cellstr(fmat{2}),idl) % region without overlap
% 
% [notInID,notInfmat, match] = compare_suj_cells(cellstr(fmat{2}),[matchr;matchl]) % region without overlap
% clear fmat
% %%fmat{1} = [matchr;notInID;matchl]
% [~ , name ] = get_parent_path([matchr;notInID;matchl],1)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

reseau1 = load(which('reseau1.mat'))
reseau1 = reseau1.reseau1

network1 = load(which('network1.mat'))
network1 = network1.network1
[a,b,c] = compare_suj_cells(reseau1,network1)

monReseau =  reseau1;%c; %reseau1
monReseau = strrep(monReseau,'.nii','_1.*mat')
fmat = get_subdir_regex_files(glm,monReseau); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reseau R&O and contrôl  nw1 sans G_front_sup 
dmodel = gdir(suj,'GLM_RSnormalised')   
glm = gdir(dmodel,'GLM2_s6'); 

ORNW = load(which('ORNW.mat'))
ORNW = ORNW.ORNW

SMN = load(which('SMN.mat'))
SMN = SMN.SMN

DAN = load(which('DAN.mat'))
DAN = DAN.DAN
DAN = addsuffixtofilenames(DAN,'_mask.*mat')
SMN = addsuffixtofilenames(SMN,'_mask.*mat')
ORNW = strrep(ORNW,'.nii','_1.*mat')


mat = [DAN;ORNW;SMN]   % 8 ;85 ;7

fmat = get_subdir_regex_files(glm,mat)
[a,b,c] = compare_suj_cells(cellstr(fmat{1}),strrep(ORNW,'.nii',''))

%manques natif space
%     {'ctx_lh_G_front_sup'}
%     {'Left-LD'           }
%     {'ctx_rh_G_front_sup'}
%     {'Right-LD'          }
%     {'MNI_Left_Botz'     }
%     {'MNI_Right_Botz'    }
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fout1 = addprefixtofilenames(dmodel,'/glm6_corr_matrix_nw1.csv');
fout2 = addprefixtofilenames(dmodel,'/glm6_pval_matrix_nw1.csv');





for nbsuj = 1:length(fmat)
    sujroi = cellstr(fmat{nbsuj});
    time_course_matrix = {};
    parfor i = 1:length(sujroi)
        l = load(sujroi{i});
        time_course_matrix{i} = l.Y;
    end
    
    time_course_matrix = cell2mat(time_course_matrix);
    
    [cor_mat, pval] = corr(time_course_matrix);
    
    %figure;imagesc(cor_mat)
    [~, roiname] = get_parent_path(sujroi);
    roiname = nettoie_dir(change_file_extension(roiname,''))

    Tcor = array2table(cor_mat, 'VariableNames', roiname);
    Tpval = array2table(pval, 'VariableNames', roiname);
    
    writetable(Tcor,fout1{nbsuj})
    writetable(Tpval,fout2{nbsuj})
    
end