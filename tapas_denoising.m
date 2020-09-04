
clear all

%% Charge all data needed if e.mat inexistent

main_dir = '/home/anna.skrzatek/data/nifti_test/ben/';
cd (main_dir)

e_PARKGAME = exam(main_dir,'PARKGAME');
e_REMINARY = exam(main_dir,'REMINARY_\w{2}_');
e = e_PARKGAME + e_REMINARY; % (3:length(e_PARKGAME)); % choose specific

e.addSerie('ACTIVATION$','run_ACTIVATION',1)
e.addSerie(        'RS$','run_RS'        ,1)

e.getSerie('run').addVolume('^f\d{3}','f',3)
e.addSerie('t1mpr_S256_0_8iso_p2$','anat_T1',1)
e.getSerie('anat_T1').addVolume('^s.*p2.nii','s',1)

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
e.getSerie('tedana').addVolume('^ts_OC','ts',1);
e.getSerie('tedana').addVolume('^dn_ts','dn_ts',1);
e.getSerie('anat').addVolume('^y'  ,'y' );
e.getSerie('anat').addVolume('^p0' ,'p0' );

e.getSerie('tedana').addVolume('^wts','wts',1);
e.getSerie('tedana').addVolume('^wdn','wdn',1);
e.getSerie('run').addVolume('^wbet','wbet',1);
e.getSerie('anat').addVolume('^wp0','wp0',1);

e.gser('tedana').addVolume('^s5wts_OC.nii','s5wts',1);
e.gser('tedana').addVolume('^s5wdn_ts_OC.nii','s5wdn',1);
e.getSerie('anat').addVolume('^wp2','wp2',1);
e.getSerie('anat').addVolume('^wp3','wp3',1);

%% Charge needed data from e object

cd /home/anna.skrzatek/data/nifti_test/PRISMA_REMINARY/
load e

fvol = e.getSerie('tedana_ACTIVATION').getVolume('^wts')
%fvol.path
avol = e.getSerie('anat').getVolume('wp0')
%avol.path
ROIvol = e.getSerie('anat').getVolume('^wp[2,3]')

par.outdir = {e.getSerie('tedana_ACTIVATION').path}
outcell = {e.getSerie('tedana_ACTIVATION').path}

%% charge and transform afni displacement parameters

e.addSerie('ACTIVATION','afni','afni_ACT',1)
e.getSerie('afni_ACT').addVolume('dfile_rall','rp1D',1)

dfile = {e.getSerie('afni_ACT').getVolume('rp').path}
output_d = {e.getSerie('tedana_ACT').path}
job_rp_afni2spm(dfile,output_d)
e.getSerie('tedana_ACT').addVolume('rp_spm','rp_spm',1)

rpvol = {e.getSerie('tedana_ACT').getVolume('rp').path}

%% get all tissue masks

nvol = length(avol);
for i = 1:nvol 
    mask{i}(1,:) = {ROIvol(i).path}; mask{i}(2,:) = {ROIvol(i+nvol).path}; 
end

par.noiseROI_mask = mask;
par.noiseROI_volume = {fvol.path};
par.rp_file = rpvol;
par.volume = {fvol.path};

par.nSlice = 630;
par.physio = 0;
par.rp = 1;
par.noiseROI = 1;
par.TR = 1600;
par.outdir = outcell

job_physio_tapas(par)