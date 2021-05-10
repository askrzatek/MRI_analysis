
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %% Charge all data needed if e.mat inexistent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Init

par.run = 1;

main_dir = '/home/anna.skrzatek/data/nifti_test/';
cd (main_dir)

e_PARKGAME = exam(main_dir,'PARKGAME.*[a,c]$');
%e_REMINARY = exam(main_dir,'REMINARY_\w{2}_');

e = e_PARKGAME; %+ e_REMINARY; % (3:length(e_PARKGAME)); % choose specific

%% Add series

e.addSerie('ACTIVATION$','run_ACTIVATION',1)
e.addSerie(        'RS$','run_RS'        ,1)

%e.addSerie('t1mpr_.*p2$','anat_T1',1)
e.addSerie('t1mpr_S256_0_8iso_p2$','anat_T1',1)

%% Add volumes

e.getSerie('run_ACT').addVolume('^v.*ACTIVATION.*nii','f',3)
e.getSerie('run_RS').addVolume('^v.*RS.*nii','f',3)

e.getSerie('anat_T1').addVolume('^v.*p2.nii','s',1)

%% old versions non-xnat
% e.getSerie('run').addVolume('^f\d{3}','f',3)
% e.addSerie('t1mpr_S256_0_8iso_p2$','anat_T1',1)
% e.getSerie('anat_T1').addVolume('^s.*p2.nii','s',1)
%%
%% Add afni outputs

e.reorderSeries('path');
e.unzipVolume(par);

e.getSerie('run').addVolume('^vtde1.nii','vtde1',1)
e.getSerie('run').addVolume('^vtde2.nii','vtde2',1)
e.getSerie('run').addVolume('^vtde3.nii','vtde3',1)

% Check if you don't need to remove some subjects (missing files)
[ec_vtd, ei_vtd] = e.removeIncomplete;
e = ec_vtd;

%% Add tedana outputs (series & volumes)

%bet
e.gser('run').addVolume('^bet.*vtde1.nii.gz','bet',1);
e.gser('run').gvol('bet').removeEmpty.unzip_and_keep(par);
e.getSerie('run').addVolume('^wbet','wbet',1);

%series
e.addSerie('ACTIVATION','tedana009.*_vtd','tedana_ACTIVATION',1);
e.addSerie('RS','tedana009.*_vtd','tedana_RS',1);

% tedana volumes
e.getSerie('tedana').addVolume('^ts_OC.*nii$','ts',1);
e.getSerie('tedana').addVolume('^dn_ts.*nii$','dn_ts',1);

e.getSerie('tedana').addVolume('^wts','wts',1);
e.getSerie('tedana').addVolume('^wdn','wdn',1);

e.gser('tedana').addVolume('^s5wts_OC.nii','s5wts',1);
e.gser('tedana').addVolume('^s5wdn_ts_OC.nii','s5wdn',1);
e.gser('tedana').addVolume('^s6wts_OC.nii','s6wts',1);
e.gser('tedana').addVolume('^s6wdn_ts_OC.nii','s6wdn',1);


%% Adding anat maps
e.getSerie('anat').addVolume('^y'  ,'y' );
e.getSerie('anat').addVolume('^p0' ,'p0' );

e.getSerie('anat').addVolume('^wp0','wp0',1);

e.getSerie('anat').addVolume('^wp2','wp2',1);
e.getSerie('anat').addVolume('^wp3','wp3',1);

e.getSerie('anat').addVolume('^rwp2','rwp2',1);
e.getSerie('anat').addVolume('^rwp3','rwp3',1);

%% 
save ('e','e')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %% Charge needed data from e object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd /home/anna.skrzatek/data/nifti_test/
load e

%% other models
%fvol = e.getSerie('tedana_ACTIVATION').getVolume('^s5wdn').toJob(0);
%fvol = e.getSerie('tedana_ACTIVATION').getVolume('^wdn').toJob(0);

%% ts + TAPAS model

% % Activation
% fvol = e.getSerie('tedana_ACTIVATION').getVolume('^wts').toJob(0);

% Resting state
fvol = e.getSerie('tedana_RS').getVolume('^wts').toJob(0);

%fvol.path % for verification

%% filter images
avol = e.getSerie('anat').getVolume('wp0');
%avol.path % for verification
ROIvol = e.getSerie('anat').getVolume('^rwp[2,3]');

% Activation params
% par.outdir = {e.getSerie('tedana_ACTIVATION').path}
% outcell = e.getSerie('run_ACT').mkdir('wts');
%%outcell = {e.getSerie('tedana_ACTIVATION').path} % old version

% Resting State params
% par.outdir = {e.getSerie('tedana_RS').path}
% outcell = e.getSerie('run_RS').mkdir('wts');

% General tedana wts dir paths (Activation + Resting State)
par.outdir = {e.getSerie('tedana').path}
outcell = e.getSerie('run').mkdir('wts');

%% charge and transform afni displacement parameters

% % Activation
% e.addSerie('ACTIVATION','afni','afni_ACT',1)
% e.getSerie('afni_ACT').addVolume('dfile_rall','rp1D',1)
% 
% dfile = {e.getSerie('afni_ACT').getVolume('rp').path}
% output_d = {e.getSerie('run_ACT').path}
% job_rp_afni2spm(dfile,output_d)
% e.getSerie('run_ACT').addVolume('rp_spm','rp_spm',1)

%rpvol = {e.getSerie('run_ACT').getVolume('rp').path}

% Resting state
e.addSerie('RS','afni','afni_RS',1)

e.getSerie('afni_RS').addVolume('dfile_rall','rp1D',1)

dfile = {e.getSerie('afni_RS').getVolume('rp').path}
output_d = {e.getSerie('run_RS').path}
job_rp_afni2spm(dfile,output_d)
e.getSerie('run_RS').addVolume('rp_spm','rp_spm',1)

%% get all tissue masks

nvol = length(avol);
par.jobname = 'tapas_get_rp';
for i = 1:nvol 
    mask{i}(1,:) = {ROIvol(i).path}; mask{i}(2,:) = {ROIvol(i+nvol).path}; 
end
% % Activation
% e.getSerie('run_ACT').addJson('^rp.*txt$','rp',1)
% rp = e.getSerie('run_ACT').getJson('rp').toJob(0);
 
% Resting state
e.getSerie('run_RS').addJson('^rp.*txt$','rp',1)
rp = e.getSerie('run_RS').getJson('rp').toJob(0);

%% par structure for tapas denoise function

par.noiseROI_mask = mask;
%par.noiseROI_volume = {fvol.path};
par.noiseROI_volume = fvol;
par.rp_file = rp;
par.volume = fvol;

par.nSlice = 60;
par.physio = 0;
par.rp = 1;
par.noiseROI = 1;
par.TR = 1.600;
par.rp_threshold = 1;
%par.rp_threshold =  %plus bas mieux c'est
%par.outdir = outcell(32:62); % Activation
par.outdir = outcell(1:31); % Resting state

par.redo = 0;
par.display = 0;
par.run = 1;
% par.job_name = 'tapas_get_rp_ACT'; % Activation
par.job_name = 'tapas_get_rp_RS'; % Resting state

job_physio_tapas( par )

%% Quality measures % names to display should be changed once we use my personally named data
%
%   [quality_measures, dR] = tapas_physio_get_movement_quality_measures(R, rHead);
%
% IN
%   R       [nScans,6]  realignment parameters estimates in mm and rad
%                       as output by SPM (x,y,z,pitch,roll,yaw)
%   rHead               head radius in mm (default: 50 mm)
%

PSNames = {'P001_NB_S1' 'P002_BM_S1' 'P001_NB_S2' 'P002_BM_S2' 'P003_SM_S1' 'P003_SM_S2' 'P007_SD_S1' 'P008_JR_S1' 'P007_SD_S2' 'P008_JR_S2' 'P023_LJ_S1' 'P025_CA_S1' 'P027_OR_S1' 'P028_PC_S1' 'P025_CA_S2' 'P023_LJ_S2' 'P027_OR_S2' 'P028_PC_S2' 'P033_DD_S1' 'P033_DD_S2' 'P039_KM_S1' 'P040_RE_S1' 'P042_RS_S1' 'P043_PD_S1' 'P039_KM_S2' 'P044_CK_S1' 'P043_PD_S2' 'P046_HJ_S1' 'P044_CK_S2' 'P047_BF_S1' 'P048_SB_S1' 'P047_BF_S2' 'P048_SB_S2'};
nsupVols_tab = [];

for r = 1 : length(rp)
    R = load (rp{r});
    rHead = 50;
    protocol = string(outcell{r}(52:55));
%     if r > 29
%         init = string(outcell{r}(end-25:end-24));
%         session = string(outcell{r}(end-21:end-20));
%     else
%         init = string(outcell{r}(end-35:end-34));
%         session = string(outcell{r}(end-21:end-20));
%     end
%     lab(r) = sprintf('%s %s %s ', protocol, init, session);
%     lab(r) = string(outcell{r}(52:end-20));
    lab(r) = string(PSNames{r});
    [quality_measures, dR] = tapas_physio_get_movement_quality_measures(R, rHead);
    mFD(r) = quality_measures.meanFD;
    sdFD(r) = std(quality_measures.FD);
    rmsMvt(r) = quality_measures.rmsMovement;
    
    nb02 = 0;
    nb05 = 0;
    nb10 = 0;
    nb25 = 0;
    
    for supVol = 1 : length(quality_measures.FD)
        if quality_measures.FD(supVol) > 0.2 
            nb02 = nb02 +1;
        end
        if quality_measures.FD(supVol) > 0.5 
            nb05 = nb05 +1;
        end
        if quality_measures.FD(supVol) > 1 
            nb10 = nb10 +1;
        end
        if quality_measures.FD(supVol) > 2.5 
            nb25 = nb25 +1;
        end
    end
    
    FDthresh = categorical({'FD > 0.2', 'FD > 0.5' 'FD > 1', 'FD > 2.5'});
    nsupVols = [nb02, nb05, nb10, nb25];
    nsupVols_tab = vertcat(nsupVols_tab, [nb02, nb05, nb10, nb25]);
    
    figure
    FDt = bar(FDthresh, nsupVols);
    title(['ACT FD: ', lab(r)], 'Interpreter', 'none')
    ylim([0, 630])
    yticks([0:15:630])
%     title(['RS FD: ', lab(r)], 'Interpreter', 'none')
%     ylim([0, 300])
%     yticks([0:5:300])
    grid on
    
    
    figure
    t = string(PSNames{r});
    brutplot = plot(quality_measures.FD);
    title(['ACT FD data by scan: ', t],'Interpreter','none')
    xlim([-20,650])
%     title(['RS FD data by scan: ', t],'Interpreter','none')
%     xlim([-20,320])
    %ylim([0,9])
end
%% sort principal subject data
[sorted_lab, sorted_index] = sort(lab); %% maybe we should use sorted lab order already in the TAPAS processing ?

%% Get quality measures (number of volumes to get rid of by subject) according to FD threshold
%s_nsupVols_tab(:,1) = (nsupVols_tab(sorted_index,1));

nsupVols_tab_id = horzcat(lab.', nsupVols_tab);
snsupVols_tab = nsupVols_tab(sorted_index,:);
sorted_nsupVols_tab_id = nsupVols_tab_id(sorted_index,:);
%sorted_yFD = [sorted_nsupVols_tab_id(:,2).'; sorted_nsupVols_tab_id(:,3).'; sorted_nsupVols_tab_id(:,4).'; sorted_nsupVols_tab_id(:,5).'];

%% Creating the txt table

[r,l] = size(sorted_nsupVols_tab_id);
nfid = fopen('nb_vols_ACT_PARK_FD.txt','w'); % ACT
%nfid = fopen('nb_vols_RS_PARK_FD.txt','w'); % RS
formatSpec = '%s;%s;%s;%s;%s\n';

tabidx = {'ID', 'FD > 0.2', 'FD > 0.5' 'FD > 1', 'FD > 2.5'};
fprintf(nfid,'%s;%s;%s;%s;%s\n',tabidx{1,:});
for k= 1:(r-1)
    fprintf(nfid,formatSpec,sorted_nsupVols_tab_id{k,1:end});
end
fprintf(nfid,'%s;%s;%s;%s;%s',sorted_nsupVols_tab_id{r,1:end});

fclose(nfid);

%% data for figure

%xFD = categorical({'FD > 0.2', 'FD > 0.5' 'FD > 1', 'FD > 2.5'});
xFD = FDthresh;
yFD = [nsupVols_tab(:,1).'; nsupVols_tab(:,2).'; nsupVols_tab(:,3).'; nsupVols_tab(:,4).'];
syFD = [snsupVols_tab(:,1).'; snsupVols_tab(:,2).'; snsupVols_tab(:,3).'; snsupVols_tab(:,4).'];

%%
figure
bpFDsub = bar(xFD, syFD, 'grouped', 'FaceColor', 'flat');
for k = 1:size(syFD,2)
    bpFDsub(k).CData = k+5;
end
ylim([0,630])
yticks([0:15:630])
title('Distribution of ACT volumes per subject according to chosen FD threshold')
% ylim([0,300])
% yticks([0:5:300])
% title('Distribution of RS volumes per subject according to chosen FD threshold')
grid on
legend(sorted_lab,'Interpreter','none')

%% Quality measures display % maybe colors to add - categories per subject

sorted_mFD = mFD(sorted_index);
sorted_sdFD = sdFD(sorted_index);
sorted_rmsMvt = rmsMvt(sorted_index);

%% sorted figure by subject name

figure
p0 = barh(sorted_mFD);
hAx = gca;
hAx.XMinorGrid = 'on';
hAx.YTick = p0.XData;
hAx.YTickLabel = sorted_lab;
hAx.TickLabelInterpreter = 'none';
hT = [];
%hT = [hT text(p1(ilab).YData, p1(ilab).XData+p1(ilab).XOffset, num2str(p1(ilab).YData.','%.3f'), 'horizontalalign','left')];
% p1.DisplayName = 'Quality measures ACT : mean FD/session/patient';
p1.DisplayName = 'Quality measures RS : mean FD/session/patient';
hAx.Title.String = p0.DisplayName;
% title('Quality measures ACT : mean FD/session/patient');
title('Quality measures RS : mean FD/session/patient');

grid on

hold on

er = errorbar(sorted_mFD, p0.XData, sorted_sdFD, 'horizontal', '+');
er.XNegativeDelta = [];
hold off

%% control figure unsorted of mean FD with sd as error bar
% figure 
% 
% p1 = barh(mFD);
% hAx = gca;
% hAx.XMinorGrid = 'on';
% hAx.YTick = p1.XData;
% hAx.YTickLabel = lab;
% hAx.TickLabelInterpreter = 'none';
% hT = [];
% %hT = [hT text(p1(ilab).YData, p1(ilab).XData+p1(ilab).XOffset, num2str(p1(ilab).YData.','%.3f'), 'horizontalalign','left')];
% % p1.DisplayName = 'Quality measures ACT : mean FD/session/patient';
% p1.DisplayName = 'Quality measures RS : mean FD/session/patient';
% hAx.Title.String = p1.DisplayName;
% grid on
% 
% hold on
% 
% er = errorbar(mFD, p1.XData, sdFD, 'horizontal', '+');
% er.XNegativeDelta = [];
% hold off

%% figure for root mean square of movement from one vol to another (sorted by subject)

figure 

p2 = barh(sorted_rmsMvt);
hAx = gca;
hAx.XMinorGrid = 'on';
hAx.YTick = p2.XData;
hAx.YTickLabel = lab;
hAx.TickLabelInterpreter = 'none';
hT = [];
%hT = [hT text(p1(ilab).YData, p1(ilab).XData+p1(ilab).XOffset, num2str(p1(ilab).YData.','%.3f'), 'horizontalalign','left')];
% p2.DisplayName = 'Quality measures ACT : total movement rms/session/patient';
p2.DisplayName = 'Quality measures RS : total movement rms/session/patient';
hAx.Title.String = p2.DisplayName;
grid on

%% number of regressors per subject, per session % NOT REALLY RELEVANT for QC but for automatic model-creation - yes man

e.addSerie('smodel_ts_tapas','model_ACT',1)
e.gser('model_ACT').addVolume('beta','beta')
e.gser('run_ACT').addVolume('wts.*.txt','multireg',1)
regcount1 = [];
regcount2 = [];

for sub = 1 : length(e)
    regpath{sub} = e(sub).gser('run_ACT').gvol('multireg').path;
    regfile = load(regpath{sub});
    [nvols nreg] = size(regfile);
    regcount1(sub)= nreg;
%     regpath{sub} = e(sub).gser('model_ACT').gvol('beta').path;
%     [nreg npath] = size(regpath{sub});
%     regcount2(sub)= nreg;
end

%% sorted figure by subject name % NOT REALLY RELEVANT
sorted_regcount = regcount1(sorted_index);

figure
preg = barh(sorted_regcount);
hAx = gca;
hAx.XMinorGrid = 'on';
hAx.YTick = preg.XData;
hAx.YTickLabel = sorted_lab;
hAx.TickLabelInterpreter = 'none';
hAx.GridAlpha
hT = [];
%hT = [hT text(p1(ilab).YData, p1(ilab).XData+p1(ilab).XOffset, num2str(p1(ilab).YData.','%.3f'), 'horizontalalign','left')];
preg.DisplayName = 'Quality measures ACT : number of regressors/session/patient';
xticks([0:50:600])
% preg.DisplayName = 'Quality measures RS : number of regressors/session/patient';
% xticks([0:10:300])
hAx.Title.String = preg.DisplayName;
grid on


%close all