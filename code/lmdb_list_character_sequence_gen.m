clear all;
clc;

rng(7610);
nchar = 15;
character = ['_', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l',...
            'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', ''''];

path = '/home/hexin/Downloads/english_words_output';
word_pre = '';
label = -1;


%% load list
list_train_path = sprintf('%s/word_train_list.txt', path);
fid = fopen(list_train_path);
line = fgetl(fid);
num = 1;

while ischar(line)
    tmp = regexp(line, '/', 'split');
    word = tmp{end-1};
    if strcmp(lower(word), word_pre)
        list{num} = sprintf('%s/%s', word, tmp{end});
    else
        word_pre = word;
        label = label + 1;
        list{num} = sprintf('%s/%s', word, tmp{end});
    end
    
    for i = 1:length(word)
        idx_label = find(character==lower(word(i)));
        list_label{num, i} = idx_label-1;
    end
    if i~=nchar
        for n = i+1:nchar
            list_label{num, n} = 0;
        end
    end
    num = num + 1;
    line = fgetl(fid);
end
fclose(fid);

%%
idx = randperm(size(list_label, 1));
num_valid = 1000;
% list_valid = list{idx(1:num_valid)};
% list_valid_label = list_label{idx(1:num_valid)};
% list_train = list{idx((num_valid+1):end)}
% list_train_label = list_label{idx((num_valid+1):end)};


%% generate validation set
fid = fopen('/home/hexin/Downloads/english_words_output//list_character_valid_lmdb.txt', 'w+');
fid_idx = fopen('/home/hexin/Downloads/english_words_output/list_character_sequence_valid_lmdb.txt', 'w+');
fprintf(fid_idx, '%d\n', num_valid);
for i = 1:num_valid
    fprintf(fid, '%s %d\n', list{idx(i)}, i-1);
    for j = 1:nchar
        fprintf(fid_idx, '%d ', list_label{idx(i), j});
    end
    fprintf(fid_idx, '\n');
end
fclose(fid);
fclose(fid_idx);




%% generate train set

fid = fopen('/home/hexin/Downloads/english_words_output/list_character_train_lmdb.txt', 'w+');
fid_idx = fopen('/home/hexin/Downloads/english_words_output/list_character_sequence_train_lmdb.txt', 'w+');
fprintf(fid_idx, '%d\n', (length(list) - num_valid));
for i = 1:(length(list) - num_valid)
    fprintf(fid, '%s %d\n', list{idx(i+num_valid)}, i-1);
    for j = 1:nchar
        fprintf(fid_idx, '%d ', list_label{idx(i+num_valid), j});
    end
    fprintf(fid_idx, '\n');
end

fclose(fid);
fclose(fid_idx);
% clear list;
% clear list_train;



%% valid set list
% list_train_path = sprintf('%s/list_word_valid.txt', path);
% fid = fopen(list_train_path)
% line = fgetl(fid);
% num = 1;
% label = -1;
% 
% while ischar(line)
%     tmp = regexp(line, '\\', 'split');
%     word = tmp{end-1};
%     if strcmp(word, word_pre)
%         list{num} = sprintf('%s/%s', word, tmp{end});
%     else
%         word_pre = word;
%         label = label + 1;
%         list{num} = sprintf('%s/%s', word, tmp{end});
%     end
%     for i = 1:length(word)
%         idx_label = find(character==word(i));
%         list_train{num, i} = idx_label;
%     end
%     if i~=nchar
%         for n = i+1:nchar
%             list_train{num, n} = 0;
%         end
%     end
%     num = num + 1;
%     line = fgetl(fid);
% end
% fclose(fid);
% 
% idx = randperm(size(list_train, 1));
% fid = fopen('../output/list_character_valid_lmdb.txt', 'w+');
% fid_idx = fopen('../output/list_character_sequence_valid_lmdb.txt', 'w+');
% fprintf(fid_idx, '%d\n', length(idx));
% for i = 1:length(idx)
%     fprintf(fid, '%s %d\n', list{idx(i)}, i-1);
%     for j = 1:nchar
%         fprintf(fid_idx, '%d ', list_train{idx(i), j});
%     end
%     fprintf(fid_idx, '\n');
% end
% 
% fclose(fid);