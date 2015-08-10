function output = data_augment(img)
        count = 1;
    %% image rotate
        angle = [-15, -10, -5, 5, 10, 15];
        n = 15;
        img = im2double(img);
        [h, w] = size(img);
        for i = 1:length(angle)
            temp = ones(h+2*n, w+2*n);
            temp(n+1:h+n, n+1:w+n) = img;
            temp = imrotate(temp, angle(i), 'bicubic');
            output{count} = temp(1.5*n:end-1.5*n, 1.5*n:end-1.5*n);
            
            count = count + 1;
        end
        
        
    %% Gaussian filter
        sigma = 20;
        gausFilter = fspecial('gaussian',  [3 3], sigma);
        temp = imfilter(img, gausFilter, 'same');


end