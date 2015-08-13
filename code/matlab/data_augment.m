function res = data_augment(img_path, save_path, param)
    input = imread(img_path);
    if size(input, 3) == 3
        input = rgb2gray(input);
    end
    if ~isa(input,'logical')
        input = input > 128;
    end
    [h, w] = size(input);
    n = 20;
    temp = ones(h+2*n, w+2*n);
    temp(n+1:h+n, n+1:w+n) = input;
    
    [height, width] = size(temp);
    input = imresize(temp, [max(height, width), max(height, width)]);
    [heigth, width] = size(input);
            
    count = 1;
    
    for iter = 1:length(param)
        wscale = param(iter).wscale;
        hscale = param(iter).hscale;
        
        w = floor(width * wscale);
        h = floor(height * hscale);
        
        temp = imresize(input, [h, w]);
        img = ones(max(w, h));
        
        if w > h
            step = floor((w - h) / 2);
            img((step+1):(step+h), :) = temp;
        else
            step = floor((h - w) / 2);
            img(:, (step+1):(step+w)) = temp;
        end
        img = imresize(img, [height, width]);
        output = im_rotate(img);
        
        for j = 1:length(output)
            tmp = regexp(img_path, '/', 'split');
            img_name = regexp(tmp{end}, '\.*', 'split');
            img_save_path = sprintf('%s/%s_%03d.png', save_path, img_name{1}, count);
            res{count} = output{j};
            output{j} = word_alignment(output{j}, 0);
            count = count + 1;           
%             imwrite(output{j}, img_save_path);
        end
    end
    
    img_save_path = sprintf('%s/%s.gif', save_path, img_name{1});
        
    for i = 1:length(res)
        imshow(res{i});
        frame = getframe(gcf);
        im=frame2im(frame);
        [I,map]=rgb2ind(im,256);
        k=i-0;
        if k==1;
            imwrite(I,map, img_save_path,'gif','Loopcount',inf,...
                'DelayTime',0.1);
        else
            imwrite(I,map,img_save_path,'gif','WriteMode','append',...
                'DelayTime',0.1);
        end
    end

end


function output = im_rotate(img)
    count = 1;
%% image rotate
    angle = [-12, -9, -6, -3, 3, 6, 9, 12];
%     angle = -12:3:12;
    n = 25;
    img = im2double(img);
    [h, w] = size(img);
    for i = 1:length(angle)
        temp = ones(h+2*n, w+2*n);
        temp(n+1:h+n, n+1:w+n) = img;
        temp = imrotate(temp, angle(i), 'bicubic');
        output{count} = temp(1.8*n:end-1.8*n, 1.8*n:end-1.8*n);
        count = count + 1;
    end
    output{count} = img;
end