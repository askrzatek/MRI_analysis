clear
clc

load e.mat
e.explore
model_name = 'model_tedana';
e.addModel(model_name,model_name);

e.getSerie('anat').addVolume('^wp0s.*','wp0')
[ec, ei] = e.removeIncomplete;

anat = ec.getSerie('anat').getVolume('wp0').getPath;
model = ec.getModel(cellstr2regex({'model_tedana'})).getPath';

%anat = repmat(anat,size(model));

coord = {
    [ -1  -9  58]  'SMA_L'
    [-36 -25  57]   'M1_L'
    [-17 -23   6]  'VIM_L'
    [  3 -68 -12] 'CB_V_R'
    [  0   0   0] 'centre'
    };

contrast = {'Real_RIGHT - Rest', 'Imaginary_RIGHT - Rest'                       % T-contrasts
    };

print_spm_figure( ...
    'model',model, ...
    'anat',anat, ...
    'coord',coord ,...
    'contrast',contrast,...
    'correction','none',...
    'threshold',0.001,...
    'nvoxel', 6,...
    'outdir',fullfile(pwd,'png'))

close all