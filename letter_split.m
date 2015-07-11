clear all;
clc;

path = '../word_train';
word_list = dir(path);
word_num = length(word_list);

% split parameter
param.pixel = 2;

% each word
for i = 4:word_num
    word_path = sprintf('%s/%s', path, word_list(i).name);
    word = word_list(i).name;
    param.nletter = length(word);
    img_list = dir(word_path);
    
    % each image 
    for i = 3:length(img_list)
        img_path = sprintf('%s/%s', word_path, img_list(i).name);
        img = imread(img_path);
        [m, n] = size(img);
        tmp = ones(1, n);
        for i = 1: m
            tmp = tmp & img(i, n);
            idx = find(tmp==0);
            letter = letter_seg(img, idx, param);
        end
    end
end