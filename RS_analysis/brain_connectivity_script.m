
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

cout = setfield(cout,'strength',strength)
cout = setfield(cout,'modularity',modularity)
cout = setfield(cout,'Elocal',Elocal)
cout = setfield(cout,'clustering',clustering)




% Stats

statGraphe = struct;
statGraphe = setfield(statGraphe,'Voi', cout.voiname')
[p, observeddifference, effectsize] = permutationTest(Eglobal(1:24), Eglobal(25:end), 200);
[psw, observeddifference, effectsize] = permutationTest(sw(1:24), sw(25:end), 200);

pV = do_statistique_test(Eglobal(1:24), Eglobal(25:end))
fields = fieldnames(cout)
fields_data = fields(6:end)

for j =1:length(fields_data)
    myfield = fields_data{j};
    
    data = getfield(cout,myfield)
        for nbr =1 :length(data) % pour chaque noeud (voi)
           pp = do_statistique_test(data(nbr,1:24), data(nbr,25:end))
           pvalue(nbr)=pp.p

        if    isnan(pvalue(nbr)) |(data(nbr,25:end) == 0  | data(nbr,1:24) == 0)
             pval(nbr) =2
        else
            pval(nbr)= permutationTest(data(nbr,1:24), data(nbr,25:end), 200); % Patients de 1:24 et contrôles 25:48
        end
        end
        %disp([pvalue' pval'])
        %pause
        statGraphe = setfield(statGraphe,myfield, pval')

        
end
 
mesIndex = (zeros(length(statGraphe.Voi),1))
myfields = fieldnames(statGraphe)

for j=2:length(myfields)

    pvals = getfield(statGraphe,myfields{j}); %#ok<*GFLD>
     %dd = str2num(cell2mat(datapemu));
    indexpemu = (pvals < 0.05)
    mesIndex = mesIndex + indexpemu

end

sout = struct;
sout = setfield(sout,'Roi',statGraphe.Voi(mesIndex>0))

for n =2:length(myfields)

    datap = getfield(statGraphe,myfields{n}); %#ok<*GFLD>
    ddata = datap(mesIndex>0)
    %ddata(mesIndex  == 0) = 1
    sout = setfield(sout,[myfields{n}],ddata)
    
end

% write_result_to_csv(sout,'resultat_reseau1_glm6.csv')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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