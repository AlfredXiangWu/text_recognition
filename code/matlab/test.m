clear all;
clc;



fox = imread('../data/20150728212540_001_fin_002_017.jpg')<128;

% Build model & skeletonize text
[fox_m, fox_rp]= autoPsm(fox);
root = find([fox_m.parent]==0);
figure
imshow(fox);
hold on;
plot([fox_m.absx]+fox_rp(1)-fox_m(root).x,[fox_m.absy]+fox_rp(2)-fox_m(root).y,'r*');




x = [fox_m.absx]+fox_rp(1)-fox_m(root).x;
y = [fox_m.absy]+fox_rp(2)-fox_m(root).y;

% [h, w] = size(fox);
% idx_l = find(x<0.5*w);
% idx_r = find(x>0.5*w);
% idx_t = find(y<0.5*h);
% idx_d = find(y>0.5*h);
% 
% figure, imshow(fox);
% hold on;
% plot(x(idx_l)-2, y(idx_l), 'b*');
% plot(x(idx_r)+2, y(idx_r), 'b*');

