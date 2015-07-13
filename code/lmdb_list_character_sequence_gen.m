clear all;
clc;

rng(7610);
nchar = 9;
character = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l',...
            'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z'];

path = '../output';
word_pre = '';
label = -1;


%% train set list
list_train_path = sprintf('%s/list_word_train.txt', path);
fid = fopen(list_train_path);
line = fgetl(fid);
num = 1;

while ischar(line)
    tmp = regexp(line, '\\', 'split');
    word = tmp{end-1};
    if strcmp(word, word_pre)
        list{num} = sprintf('%s/%s %d\n', word, tmp{end}, label);
    else
        word_pre = word;
        label = label + 1;
        list{num} = sprintf('%s/%s %d\n', word, tmp{end}, label);
    end
    for i = 1:length(word)
        idx_label = find(character==word(i));
        list_train{num, i} = idx_label;
    end
    if i~=nchar
        for n = i+1:nchar
            list_train{num, n} = 0;
        end
    end
    num = num + 1;
    line = fgetl(fid);
end
fclose(fid);

idx = randperm(size(list_train, 1));
fid = fopen('../output/list_character_train_lmdb.txt', 'w+');
fid_idx = fopen('../output/list_character_sequence_train_lmdb.txt', 'w+');
for i = 1:length(idx)
    fprintf(fid, '%s', list{idx(i)});
    for j = 1:nchar
        fprintf(fid_idx, '%d ', list_train{idx(i), j});
    end
    fprintf(fid_idx, '\n');
end

fclose(fid);
clear list;
clear list_train;



%% valid set list
list_train_path = sprintf('%s/list_word_valid.txt', path);
fid = fopen(list_train_path);
line = fgetl(fid);
num = 1;
label = -1;

while ischar(line)
    tmp = regexp(line, '\\', 'split');
    word = tmp{end-1};
    if strcmp(word, word_pre)
        list{num} = sprintf('%s/%s %d\n', word, tmp{end}, label);
    else
        word_pre = word;
        label = label + 1;
        list{num} = sprintf('%s/%s %d\n', word, tmp{end}, label);
    end
    for i = 1:length(word)
        idx_label = find(character==word(i));
        list_train{num, i} = idx_label;
    end
    if i~=nchar
        for n = i+1:nchar
            list_train{num, n} = 0;
        end
    end
    num = num + 1;
    line = fgetl(fid);
end
fclose(fid);

idx = randperm(size(list_train, 1));
fid = fopen('../output/list_character_valid_lmdb.txt', 'w+');
fid_idx = fopen('../output/list_character_sequence_valid_lmdb.txt', 'w+');
for i = 1:length(idx)
    fprintf(fid, '%s', list{idx(i)});
    for j = 1:nchar
        fprintf(fid_idx, '%d ', list_train{idx(i), j});
    end
    fprintf(fid_idx, '\n');
end

fclose(fid);