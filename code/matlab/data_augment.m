function res = data_augment(img_path, save_path, param)
    input = imread(img_path) > 128;
    [h, w] = size(input);
    n = 15;
    temp = ones(h+2*n, w+2*n);
    temp(n+1:h+n, n+1:w+n) = input;
    img = im2bw(temp);
    [height, width] = size(img);
    img = ~img;
    
    
    [psm, rp]= autoPsm(img);
    root = find([psm.parent]==0);

    x = [psm.absx]+rp(1)-psm(root).x;
    y = [psm.absy]+rp(2)-psm(root).y;
    r = psm(1).r*2;
    
    count = 1;
    
    for i = 1: length(param)
        [x_trans, y_trans] = point_transform(x, y, height, width, param(i));
        img_trans = gen_img(img, x, y, x_trans, y_trans, r);
        img_trans = ~img_trans;
        output = im_rotate(img_trans);
        
        for j = 1:length(output)
            tmp = regexp(img_path, '/', 'split');
            img_name = regexp(tmp{end}, '\.*', 'split');
            img_save_path = sprintf('%s/%s_%03d.png', save_path, img_name{1}, count);
            count = count + 1;
            imwrite(output{j}, img_save_path);
        end
    end     
    
    res = 1;
end


function output = im_rotate(img)
    count = 1;
%% image rotate
    angle = [-15, -10, -5, 5, 10, 15];
    n = 20;
    img = im2double(img);
    [h, w] = size(img);
    for i = 1:length(angle)
        temp = ones(h+2*n, w+2*n);
        temp(n+1:h+n, n+1:w+n) = img;
        temp = imrotate(temp, angle(i), 'bicubic');
        output{count} = temp(2*n:end-2*n, 2*n:end-2*n);
        count = count + 1;
    end
    output{count} = img;
end