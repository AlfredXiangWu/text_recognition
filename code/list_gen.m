clear all;
clc;

path = '/home/hexin/Downloads/english_words_output/word_train';
word_list = dir('/home/hexin/Downloads/english_words_output/word_train');

list = fopen('/home/hexin/Downloads/english_words_output/word_train_list.txt', 'w+');


for nword = 3:length(word_list)
    word_path = sprintf('%s/%s', path, word_list(nword).name);
    word = word_list(nword).name;
    param.nletter = length(word);
    img_list = dir(word_path);
    
    % save file path
    output_path = sprintf('/home/hexin/Downloads/english_words_output/word_train/%s', word);
    
    % each image 
    for i = 3:length(img_list)
        img_name = img_list(i).name;
        img_path = sprintf('%s/%s', word_path, img_name);
        fprintf(list, '%s\n', img_path);
    end
end