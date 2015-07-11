clear all;
clc;

path = '../word_train';
error_log = fopen('../output/error_log.txt', 'w+');
word_list = dir(path);
word_num = length(word_list);

% split parameter
letter_flag = 0;
param.pixel = 2;

% each word
for i = 4:word_num
    word_path = sprintf('%s/%s', path, word_list(i).name);
    word = word_list(i).name;
    param.nletter = length(word);
    img_list = dir(word_path);
    
    % save file path
    output_path = sprintf('../output/word_train/%s', word);
    mkdir(output_path);
    
    % each image 
    for i = 3:length(img_list)
        img_name = img_list(i).name;
        if strcmp(img_name(end-2:end), 'txt')
            continue;
        end     
        img_path = sprintf('%s/%s', word_path, img_name);
        
        if ~exist(img_path, 'file')
            fprintf(error_log, '%s is not exsit\n', img_path);
            continue;
        end
        
        img = imread(img_path);
        [m, n] = size(img);
        tmp = ones(1, n);
        for i = 1: m
            tmp = tmp & img(i, :);
        end
        idx = find(tmp==0);
        if letter_flag
            letter = letter_seg(img, idx, param);
        else
             save_path = sprintf('%s/%s', output_path, img_name);
             y_start = max(1, idx(1) - param.pixel);
             y_end = min(n, idx(end) + param.pixel);
             save_img = img(:, y_start:y_end)*255;
             fprintf('%s is done!\n', save_path);
             imwrite(save_img, save_path);
             clear save_img;
        end
    end
end

fclose(error_log);