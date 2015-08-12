clear all;
clc;

% original = imread('../../data/20150728212540_001_fin_002_017.jpg')<128;
original = rgb2gray(imread('./PSM/fox.png'))<128;
[h, w] = size(original);
n = 15;
temp = zeros(h+2*n, w+2*n);
temp(n+1:h+n, n+1:w+n) = original;
img = temp >0.5;
[height, width] = size(img);

%% get stoken
[psm, rp]= autoPsm(img);
root = find([psm.parent]==0);
figure
imshow(img);

x = [psm.absx]+rp(1)-psm(root).x;
y = [psm.absy]+rp(2)-psm(root).y;

%% draw 
% for i = 1:length(fox_m)
%      if isempty(fox_m(i).children)
% %     if isempty(fox_m(i).children) || (fox_m(i).parent == 0)
%         continue;
%     end
%     
%     plot([x(i) x(fox_m(i).children)], [y(i) y(fox_m(i).children)], 'b-', 'linewidth', 5);
% end


%% transform
% 1
param.top = -5;
param.bottom = 3;
param.left = 3;
param.right = 3;
r = psm(1).r*2;

[x_trans, y_trans] = point_transform(x, y, height, width, param);
output = gen_img(img, x, y, x_trans, y_trans, r);
output = ~output;
imwrite(output, 'fox1.png');

% 2
param.top = -2;
param.bottom = -2;
param.left = 2;
param.right = 3;
r = psm(1).r*2;

[x_trans, y_trans] = point_transform(x, y, height, width, param);
output = gen_img(img, x, y, x_trans, y_trans, r);
output = ~output;
imwrite(output, 'fox2.png');

% 3
param.top = 4;
param.bottom = -4;
param.left = 3;
param.right = -3;
r = psm(1).r*2;

[x_trans, y_trans] = point_transform(x, y, height, width, param);
output = gen_img(img, x, y, x_trans, y_trans, r);
output = ~output;
imwrite(output, 'fox3.png');

% 4
param.top = 4;
param.bottom = -4;
param.left = -4;
param.right = 4;
r = psm(1).r*2;

[x_trans, y_trans] = point_transform(x, y, height, width, param);
output = gen_img(img, x, y, x_trans, y_trans, r);
output = ~output;
imwrite(output, 'fox4.png');



