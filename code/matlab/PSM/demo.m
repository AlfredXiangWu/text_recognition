% Demo of Inkblot Part-Structured Models
%
% Reference code for this paper:
% N. Howe, Part-Structured Inkball Models for One-Shot Handwritten Word Spotting.
% International Conference on Document Analysis and Recognition, August 2013.

% Read in images
fox = rgb2gray(imread('fox.png'))<128;
text = rgb2gray(imread('pangrams.png'))<128;

% Build model & skeletonize text
[fox_m, fox_rp]= autoPsm(fox);
root = find([fox_m.parent]==0);
sktext = bwmorph(text,'thin',inf);
figure
imshow(fox);
hold on;
plot([fox_m.absx]+fox_rp(1)-fox_m(root).x,[fox_m.absy]+fox_rp(2)-fox_m(root).y,'r*');

% Fit model
try
    % attempt gpu version if available
    [dtsq,loc] = psmFit_gpu(fox_m,sktext,[8 8]);  % Third argument should be greater than maximum node displacement
catch
    % else fall back to cpu code
    [dtsq,loc] = psmFit(fox_m,sktext);
end;

% display best location
[~,y,x] = min2d(dtsq);
figure
imshow(text)
hold on
plot(loc{y,x}(1,:),loc{y,x}(2,:),'r*')
