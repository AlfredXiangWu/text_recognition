clear all;
clc;
% original = rgb2gray(imread('./PSM/fox.png'))<128;
% [h, w] = size(original);
% n = 15;
% temp = zeros(h+2*n, w+2*n);
% temp(n+1:h+n, n+1:w+n) = original;
% img = temp >0.5;
% [height, width] = size(img);

% %% get stoken
% [psm, rp]= autoPsm(img);
% root = find([psm.parent]==0);
% figure
% imshow(img);
% 
% x = [psm.absx]+rp(1)-psm(root).x;
% y = [psm.absy]+rp(2)-psm(root).y;

%% draw 
% for i = 1:length(fox_m)
%      if isempty(fox_m(i).children)
% %     if isempty(fox_m(i).children) || (fox_m(i).parent == 0)
%         continue;
%     end
%     
%     plot([x(i) x(fox_m(i).children)], [y(i) y(fox_m(i).children)], 'b-', 'linewidth', 5);
% end


%% configuration
% 1
param(1).top = -5;
param(1).bottom = 5;
param(1).left = -3;
param(1).right = 3;

% 2
param(2).top = -5;
param(2).bottom = 5;
param(2).left = 3;
param(2).right = -3;

% 3
param(3).top = 4;
param(3).bottom = -4;
param(3).left = 3;
param(3).right = -3;

% 4
param(4).top = 4;
param(4).bottom = -4;
param(4).left = -4;
param(4).right = 4;

% 5
param(4).top = 0;
param(4).bottom = -3;
param(4).left = 0;
param(4).right = -3;

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
