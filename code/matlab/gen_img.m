function output = gen_img(img, x, y, x_trans, y_trans, r)
    [height, width] = size(img);
    output = zeros(height, width);
    for i = 1:length(x)
        xstart = max(x(i) - r, 1);
        xend = min(x(i) + r, width);
        ystart = max(y(i) - r, 1);
        yend = min(y(i) + r, height);
        temp = img(ystart:yend, xstart:xend);

        xstart = max(x_trans(i) - r, 1);
        xend = min(x_trans(i) + r, width);
        ystart = max(y_trans(i) - r, 1);
        yend = min(y_trans(i) + r, height);
        output(ystart:yend, xstart:xend) = temp;
    end

    se = strel('square', r);
    bw = imdilate(output, se);
%     se = strel('square', r);
    bw = imerode(bw, se);
    bw = imerode(bw, se);
%     bw = imdilate(bw, se);
    output = bw;
end