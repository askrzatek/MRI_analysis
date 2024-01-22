clear all 
clc
suj = getSubjects;
[patients, controls]= separate(suj)


suj = [patients;controls];



pool = cell(length(suj),1)
pool(1:length(patients),1) = {'1'}
pool(length(patients)+1:length(suj),1)  = {'2'}
[~,sujname] = get_parent_path(suj,1);
    
dmodel = gdir(suj,'GLM_Native')%'^GLM_RSnormalised')  %

corrMatrix = gfile(dmodel,'^glm6_corr_matrix_nw1.*csv')

Xdan = 1:8
Xonw = 9:56
Xrnw =21:89
Xsmn = 90:96

Xornw = 9:89
for i = 1:length(corrMatrix)
    
    tablemtx = readtable(corrMatrix{i});
    voiname = tablemtx.Properties.VariableNames
    mtx = table2array(tablemtx);
    a = threshold_absolute(mtx,-1);
    z = atanh(a)
    w = abs(z)
    
    dan = w(Xdan,Xdan)
    onw = w(Xonw,Xonw)
    rnw = w(Xrnw,Xrnw)
    smn = w(Xsmn,Xsmn)
    ornw = w(Xornw,Xornw)
    
    zdan(i,1) = mean_connection_strength(dan)
    zonw(i,1) = mean_connection_strength(onw)
    zrnw(i,1) = mean_connection_strength(rnw)
    zsmn(i,1) = mean_connection_strength(smn)
    
    zornw(i,1) = mean_connection_strength(ornw)
    
    % a = sum(sum(dan))/((8*(7)/2))
    % W = weight_conversion(w, 'binarize');
    
    dan_smn = w(Xdan,Xsmn)
    onw_dan = w(Xonw,Xdan)
    rnw_dan = w(Xrnw,Xdan)
    
    onw_smn = w(Xonw,Xsmn)
    rnw_smn = w(Xrnw,Xsmn)
    onw_rnw = w(Xonw,Xrnw)
    
    ornw_smn = w(Xornw,Xsmn)
    ornw_dan = w(Xornw,Xdan)
    
    
    zdan_smn(i,1) =  mean_strength_inter(dan_smn)
    zonw_dan(i,1) =  mean_strength_inter(onw_dan)
    zrnw_dan(i,1) =  mean_strength_inter(rnw_dan)
    zonw_smn(i,1) =  mean_strength_inter(onw_smn)
    zrnw_smn(i,1) =  mean_strength_inter(rnw_smn)
    
    zornw_smn(i,1) =  mean_strength_inter(ornw_smn)
    zornw_dan(i,1) =  mean_strength_inter(ornw_dan)


    % sw(i) = small_worldness(W,100)    
    
    Integrity(:,i) = (strengths_und(w))'; % C'est la somme des poids des liens dans un noeud 
    
    
%     modularity(:,i) = modularity_und(W); % c'est le nombre de groupe ou de communauté
%     Eglobal(i) = efficiency_bin(W);   
%     Elocal(:,i)=  efficiency_bin(W,1);
%     clustering(:,i) = clustering_coef_bu(W);
    

end



cout = struct;
cout = setfield(cout,'suj',sujname)
cout = setfield(cout,'pool',pool)
%cout = setfield(cout,'voiname',voiname)
cout = setfield(cout,'Z_DAN',zdan)
cout = setfield(cout,'Z_ONW',zonw)
cout = setfield(cout,'Z_RNW',zrnw)
cout = setfield(cout,'Z_SMN',zsmn)
cout = setfield(cout,'Z_ORNW',zornw)


cout = setfield(cout,'Z_DAN_SMN',zdan_smn)
cout = setfield(cout,'Z_ONW_DAN',zonw_dan)
cout = setfield(cout,'Z_RNW_DAN',zrnw_dan)

cout = setfield(cout,'Z_ONW_SMN', zonw_smn)
cout = setfield(cout,'Z_RNW_SMN', zrnw_smn)
cout = setfield(cout,'Z_ORNW_SMN',zornw_smn )
cout = setfield(cout,'Z_ORNW_DAN',zornw_dan )

write_result_to_csv(cout,'inter_nw_mesures_natif.csv')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



fields = fieldnames(cout)
fields_data = fields(3:end)

for j =1:length(fields_data)
    myfield = fields_data{j};
    
    data = getfield(cout,myfield)
        
    pp = do_statistique_test(data(1:length(patients)), data(length(patients)+1:end));
    pvalue(j)=pp.p

     
    pval(j)= permutationTest(data(1:length(patients)), data(length(patients)+1:end), 200) % Patients de 1:24 et contrôles 25:48
       
        %pause
       % statGraphe = setfield(statGraphe,myfield, pval')    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for nbr =1:(96)
    
    
    pp = do_statistique_test(Integrity(nbr,1:length(patients)), Integrity(nbr,length(patients)+1:end))
    zpvalue(nbr)=pp.p

       
    zpval(nbr)= permutationTest(Integrity(nbr,1:length(patients)), Integrity(nbr,length(patients)+1:end), 200); % Patients de 1:24 et contrôles 25:48
      
        
end

name = strrep(voiname,'VOI_VOI_w','')
name = strrep(voiname,'VOI_VOI_','')
name = strrep(name,'VOI_w','')
name = strrep(name,'_1','')

tab = array2table(Integrity','VariableNames',name')
po = table(pool)
a = zpval <0.05
b = zpvalue <0.05
an = voiname(a')      %% {'VOI_rcLeft_VIIIa_1'}    {'VOI_rcRight_VIIIa_1'}
bn = voiname(b')

tab = [po tab]
writetable(tab,'integrity_nw1_natif.csv')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%


resultat = readtable('resultat_reseau1_glm6.csv');
voi = resultat.Roi
name = strrep(voi,'VOI_ra','')
name = strrep(name,'VOI_rs','')
name = strrep(name,'VOI_rt','')
name = strrep(name,'VOI_rh','')
name = strrep(name,'VOI_rc','')
name = strrep(name,'VOI_rw','')
name = strrep(name,'VOI_ry','')
name = strrep(name,'_1','')


metric_network = resultat.Properties.VariableNames
metric_network = metric_network(2:end-1)


for nbrV =  1:length(metric_network)
    data_metric = getfield(resultat,metric_network{nbrV})
    mesVoi = voi(data_metric<0.05)
    mesnoms =name((data_metric<0.05))
    
    metric = struct
    metric = setfield(metric,'suj',sujname)
    metric = setfield(metric,'pool',pool)
    for nbrvoi = 1: length(mesVoi)
        indexVoi = find(strcmp(voiname,mesVoi{nbrvoi}))
        
        if strcmp(metric_network{nbrV}, 'degree')
            metric = setfield(metric,strrep(mesnoms{nbrvoi},'-','_'),deg(indexVoi,:)')
        else
            data_metric2 = eval(metric_network{nbrV})
             metric = setfield(metric,strrep(mesnoms{nbrvoi},'-','_'),data_metric2(indexVoi,:)')
        end
    end

   write_result_to_csv(metric,['GLM6_' metric_network{nbrV} '.csv']);
end
