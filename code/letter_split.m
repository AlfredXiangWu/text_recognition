clear all;
clc;

path = '/home/hexin/Downloads/english_words';
error_log = fopen('/home/hexin/Downloads/english_words_output/train_error_log.txt', 'w+');
word_list = dir(path);
word_num = length(word_list);

output_list = fopen('/home/hexin/Downloads/english_words_output/word_train_list.txt', 'w+');


% split parameter
letter_flag = 0;
param.pixel = 2;

% each word
for nword = 4:word_num
    word_path = sprintf('%s/%s', path, word_list(nword).name);
    word = word_list(nword).name;
    
    tmp = regexp(word, ' ', 'split');
    if length(tmp)~=1
       word = tmp{1};
       for nn = 1:length(tmp)-1
          word = sprintf('%s_%s', word, tmp{nn+1}); 
       end
    end
    
    param.nletter = length(word);
    img_list = dir(word_path);
    
    % save file path
    output_path = sprintf('/home/hexin/Downloads/english_words_output/word_train/%s', word);
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
        for j = 1: m
            tmp = tmp & img(j, :);
        end
        idx = find(tmp==0);
        if isempty(idx)
            continue;
        end
            
        if letter_flag
            letter = letter_seg(img, idx, param);
        else
             save_path = sprintf('%s/%d.png', output_path, i-2);
             fprintf(output_list, '%s\n', save_path);
             
             if length(idx) > 0.5*m
                 save_img = uint8(255*ones(m, 2*m));                      
                 y_start = max(1, idx(1) - param.pixel);
                 y_end = min(n, idx(end) + param.pixel);
                 save_img(:, 1:(y_end-y_start+1)) = img(:, y_start:y_end);
             else
                 y_start = max(1, idx(1) - param.pixel);
                 y_end = min(n, idx(end) + param.pixel);
                 save_img = img(:, y_start:y_end);
             end
             save_img = imresize(save_img, [50 100]);
             fprintf('%s is done!\n', save_path);
             imwrite(save_img, save_path);
             clear save_img;
        end
    end
end

fclose(error_log);