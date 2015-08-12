function [x_trans, y_trans] = point_transform(x, y, height, width, param)
    top = param.top;
    bottom = param.bottom;
    left = param.left;
    right = param.right;
    
    x_trans = x;
    y_trans = y;
    
    % top-bottom
    idx_t = find(y < 0.6*height);
    idx_b = find(y > 0.4*height);
    
    y_trans(idx_t) = max(y_trans(idx_t) +top, 1);
    y_trans(idx_b) = min(y_trans(idx_b) + bottom, height);
    
    % left-right
    idx_l = find(x < 0.6*width);
    idx_r = find(x > 0.4*width);
    
    x_trans(idx_l) = max(x_trans(idx_l) + left, 1);
    x_trans(idx_r) = min(x_trans(idx_r)+ right, width);

end