clear all 
clc
suj = getSubjects;
[patients, controls]= separate(suj)


suj = [patients;controls];



pool = cell(length(suj),1)
pool(1:length(patients),1) = {'1'}
pool(length(patients)+1:length(suj),1)  = {'2'}
[~,sujname] = get_parent_path(suj,1);
    
dmodel = gdir(suj,'^GLM_NativeSapce')  %

corrMatrix = gfile(dmodel,'^glm6_corr_matrix_nw1.*csv')



for i = 1:length(corrMatrix)
    
    tablemtx = readtable(corrMatrix{i});
    voiname = tablemtx.Properties.VariableNames
    mtx = table2array(tablemtx);
    w = abs(mtx)
    w = threshold_absolute(w,0);
    W = ECOfilter(w,0)
    wW = w
    wW(W==0)=0
    % W = weight_conversion(w, 'binarize');
    
    
    deg(:,i) = degrees_und(W);
    sw(i) = small_worldness(W,100)    
    
    strength(:,i) = (strengths_und(wW))'; % C'est la somme des poids des liens dans un noeud 
    
    
    modularity(:,i) = modularity_und(W); % c'est le nombre de groupe ou de communauté
    Eglobal(i) = efficiency_bin(W);   
    
    Elocal(:,i)=  efficiency_bin(W,1);
    
    clustering(:,i) = clustering_coef_bu(W);

end


cout = struct;
cout = setfield(cout,'suj',sujname')
cout = setfield(cout,'pool',pool)
cout = setfield(cout,'voiname',voiname)
cout = setfield(cout,'sw',sw)
cout = setfield(cout,'Eglobal',Eglobal)
cout = setfield(cout,'degree',deg)
