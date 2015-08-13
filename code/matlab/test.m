clear all;
clc;

% 1
param(1).wscale = 1.5;
param(1).hscale = 1;

% 2
param(2).wscale = 1;
param(2).hscale = 1.5;

% 3
param(3).wscale = 2;
param(3).hscale = 1;

% 4
param(4).wscale = 1;
param(4).hscale = 2;


%% process
path = '../../data';
% img_path= '../../data/fox.png';
save_path = '../../results';

subdir = dir(path);
for i  =1:length(subdir)
    if subdir(i).isdir
        continue;
    end
    
    img_path = sprintf('%s/%s', path, subdir(i).name);
    res = data_augment(img_path, save_path, param);
end
