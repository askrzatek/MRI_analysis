
clear all

%% Charge all data needed if e.mat inexistent

main_dir = '/home/anna.skrzatek/data/nifti_test/';
cd (main_dir)

e_PARKGAME = exam(main_dir,'PARKGAME.*[a,c]$');
%e_REMINARY = exam(main_dir,'REMINARY_\w{2}_');
e = e_PARKGAME %+ e_REMINARY; % (3:length(e_PARKGAME)); % choose specific

%
e.addSerie('ACTIVATION$','run_ACTIVATION',1)
e.addSerie(        'RS$','run_RS'        ,1)
% 
e.getSerie('run_ACT').addVolume('^v.*ACTIVATION.*nii','f',3)
e.getSerie('run_RS').addVolume('^v.*RS.*nii','f',3)
%e.addSerie('t1mpr_.*p2$','anat_T1',1)
e.addSerie('t1mpr_S256_0_8iso_p2$','anat_T1',1)
e.getSerie('anat_T1').addVolume('^v.*p2.nii','s',1)

% e.getSerie('run').addVolume('^f\d{3}','f',3)
% e.addSerie('t1mpr_S256_0_8iso_p2$','anat_T1',1)
% e.getSerie('anat_T1').addVolume('^s.*p2.nii','s',1)

e.reorderSeries('path');
e.unzipVolume(par);
e.getSerie('run').addVolume('^vtde1.nii','vtde1',1)
e.getSerie('run').addVolume('^vtde2.nii','vtde2',1)
e.getSerie('run').addVolume('^vtde3.nii','vtde3',1)

[ec_vtd, ei_vtd] = e.removeIncomplete;
e = ec_vtd;

e.gser('run').addVolume('^bet.*vtde1.nii.gz','bet',1);
e.gser('run').gvol('bet').removeEmpty.unzip_and_keep(par);
e.addSerie('ACTIVATION','tedana009.*_vtd','tedana_ACTIVATION',1);
e.addSerie('RS','tedana009.*_vtd','tedana_RS',1);
e.getSerie('tedana').addVolume('^ts_OC.*nii$','ts',1);
e.getSerie('tedana').addVolume('^dn_ts.*nii$','dn_ts',1);
e.getSerie('anat').addVolume('^y'  ,'y' );
e.getSerie('anat').addVolume('^p0' ,'p0' );

e.getSerie('tedana').addVolume('^wts','wts',1);
e.getSerie('tedana').addVolume('^wdn','wdn',1);
e.getSerie('run').addVolume('^wbet','wbet',1);
e.getSerie('anat').addVolume('^wp0','wp0',1);

e.gser('tedana').addVolume('^s5wts_OC.nii','s5wts',1);
e.gser('tedana').addVolume('^s5wdn_ts_OC.nii','s5wdn',1);

e.gser('tedana').addVolume('^s6wts_OC.nii','s6wts',1);
e.gser('tedana').addVolume('^s6wdn_ts_OC.nii','s6wdn',1);


e.getSerie('anat').addVolume('^wp2','wp2',1);
e.getSerie('anat').addVolume('^wp3','wp3',1);

save ('e','e')

%% Charge needed data from e object

cd /home/anna.skrzatek/data/nifti_test/
load e

%fvol = e.getSerie('tedana_ACTIVATION').getVolume('^s5wdn').toJob(0);
%fvol = e.getSerie('tedana_ACTIVATION').getVolume('^wdn').toJob(0);
%fvol = e.getSerie('tedana_ACTIVATION').getVolume('^wts').toJob(0);
fvol = e.getSerie('tedana_RS').getVolume('^wts').toJob(0);
%fvol.path
avol = e.getSerie('anat').getVolume('wp0');
%avol.path
ROIvol = e.getSerie('anat').getVolume('^wp[2,3]');

% par.outdir = {e.getSerie('tedana_ACTIVATION').path}
% outcell = e.getSerie('run_ACT').mkdir('wts');
%outcell = {e.getSerie('tedana_ACTIVATION').path}

% par.outdir = {e.getSerie('tedana_RS').path}
% outcell = e.getSerie('run_RS').mkdir('wts');

par.outdir = {e.getSerie('tedana').path}
outcell = e.getSerie('run').mkdir('wts');

%% charge and transform afni displacement parameters

% e.addSerie('ACTIVATION','afni','afni_ACT',1)
% e.getSerie('afni_ACT').addVolume('dfile_rall','rp1D',1)
% 
% dfile = {e.getSerie('afni_ACT').getVolume('rp').path}
% output_d = {e.getSerie('run_ACT').path}
% job_rp_afni2spm(dfile,output_d)
% e.getSerie('run_ACT').addVolume('rp_spm','rp_spm',1)
% 
% %rpvol = {e.getSerie('run_ACT').getVolume('rp').path}

ei.addSerie('RS','afni','afni_RS',1)
ei.getSerie('afni_RS').addVolume('dfile_rall','rp1D',1)

dfile = {ei.getSerie('afni_RS').getVolume('rp').path}
output_d = {ei.getSerie('run_RS').path}
job_rp_afni2spm(dfile,output_d)
ei.getSerie('run_RS').addVolume('rp_spm','rp_spm',1)

%% get all tissue masks

nvol = length(avol);
par.jobname = 'tapas_get_rp';
for i = 1:nvol 
    mask{i}(1,:) = {ROIvol(i).path}; mask{i}(2,:) = {ROIvol(i+nvol).path}; 
end
% e.getSerie('run_ACT').addJson('^rp.*txt$','rp',1)
% rp = e.getSerie('run_ACT').getJson('rp').toJob(0);

ei.getSerie('run_RS').addJson('^rp.*txt$','rp',1)
rp = ei.getSerie('run_RS').getJson('rp').toJob(0);

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
%par.rp_threshold =  %plus bas mieux c'est
par.outdir = outcell

par.redo = 1;
par.display = 0;
par.run = 1;
par.job_name = 'tapas_get_rp_P046';

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
%     label(r) = sprintf('%s %s %s ', protocol, init, session);
    label(r) = string(outcell{r}(52:end-20));
    [quality_measures, dR] = tapas_physio_get_movement_quality_measures(R, rHead);
    mFD(r) = quality_measures.meanFD;
end

%% Quality measures display % maybe colors to add - categories per subject
%label
figure 

p1 = barh(mFD);
hAx = gca;
hAx.XMinorGrid = 'on';
hAx.YTick = p1.XData;
hAx.YTickLabel = label;
hT = [];
hT = [hT text(p1(ilabel).YData, p1(ilabel).XData+p1(ilabel).XOffset, num2str(p1(ilabel).YData.','%.3f'), 'horizontalalign','left')];
p1.DisplayName = 'Quality measures : mean FD/session/patient';
hAx.Title.String = p1.DisplayName;

%close all