function do_REMINARY_PARKGAME_architecture(rootdesti)
%main_dir = fullfile('/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek','nifti');

% e_PARKGAME = exam(main_dir,'PARKGAME');
% e_REMINARY = exam(main_dir,'REMINARY_\w{2}_');
e_REMINARYS1 = exam(main_dir,'REMINARY_\w{2}_.*1$');
e_REMINARYS2 = exam(main_dir,'REMINARY_\w{2}_.*2$');
e_PARKGAMES2 = exam(main_dir,'PARKGAME.*2$');
e_PARKGAMES1 = exam(main_dir,'PARKGAME.*1$');

e = {e_PARKGAMES1, e_PARKGAMES2, e_REMINARYS1, e_REMINARYS2};

% for i = 1:length(e)
%     %e{i}.explore
%     %'REMINARY_\w{2}_.*1$'
%     e{i}.addSerie('model_meica$','contrasts',1)
% 
%     e{i}.getSerie('contrasts').addVolume('^con_0008','REAL_L',1)
%     e{i}.getSerie('contrasts').addVolume('^con_0009','REAL_R',1)
%     e{i}.getSerie('contrasts').addVolume('^con_0010','IMA_L',1)
%     e{i}.getSerie('contrasts').addVolume('^con_0011','IMA_R',1)
% 
%     %e{i}.explore
% end

dirstat = r_mkdir(main_dir, 'secondlevel_ACTIVATION_test');
dirgroup = r_mkdir(char(dirstat), {'PARKGAME', 'REMINARY'});


nsess = 2;
ncond = nsess*length(dirgroup);
igroup = 0;

for i = 1:nsess:ncond
    
   s1 = e{i}.getSerie('contrasts').toJob;
   cons1 = get_subdir_regex_files(s1, '^con');
   cons1{1}(1,:)
   s2 = e{i+1}.getSerie('contrasts').toJob;
   cons2 = get_subdir_regex_files(s2, '^con');

   igroup = igroup+1;
   dirsess = r_mkdir(dirgroup{igroup}, {'S1', 'S2'});
   for k = 1 : length(cons1)
       for kk = 1 : length(cons1{k}(:,k))
           cmd='';
           cmd =  sprintf('%s\nln -f -s %s %s',cmd, cons1{k}(kk,:), dirsess{i});
           unix(cmd);
       end
   end
   for k = 1 : length(cons2)
       cmd='';
       cmd =  sprintf('%s\nln -f -s %s %s',cmd, s2{k} , dirsess{i+1});
   end
end

end
