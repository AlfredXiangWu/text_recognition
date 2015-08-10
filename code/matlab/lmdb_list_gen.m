clear all;
clc;

path = '../output';
word_pre = '';
label = -1;

%% train set list
list_train_path = sprintf('%s/list_word_train.txt', path);
fid = fopen(list_train_path);
line = fgetl(fid);
list_train_output = fopen('../output/list_word_train_lmdb.txt', 'w+');


while ischar(line)
    tmp = regexp(line, '\\', 'split');
    word = tmp{end-1};
    if strcmp(word, word_pre)
        fprintf(list_train_output, '%s/%s %d\n', word, tmp{end}, label);
    else
        word_pre = word;
        label = label + 1;
        fprintf(list_train_output, '%s/%s %d\n', word, tmp{end}, label);
    end
    line = fgetl(fid);
end

fclose(fid);
fclose(list_train_output);

%% valid set list
word_pre = '';
label = -1;
list_valid_path = sprintf('%s/list_word_valid.txt', path);
fid = fopen(list_valid_path);
line = fgetl(fid);
list_valid_output = fopen('../output/list_word_valid_lmdb.txt', 'w+');


while ischar(line)
    tmp = regexp(line, '\\', 'split');
    word = tmp{end-1};
    if strcmp(word, word_pre)
        fprintf(list_valid_output, '%s/%s %d\n', word, tmp{end}, label);
    else
        word_pre = word;
        label = label + 1;
        fprintf(list_valid_output, '%s/%s %d\n', word, tmp{end}, label);
    end
    line = fgetl(fid);
end

fclose(fid);
fclose(list_valid_output);