function [res] = word_alignment(img, step)
    [height, width] = size(img);
    
    binary = img > 0.8;
    
%     idx = [];
%     v_avg = mean(binary, 1);   
%     idx = find(v_avg < 0.98);
%     if isempty(idx)
%         res = binary;
%         position.hstart = 0;
%         position.hend = height;
%         position.wstart = 0;
%         position.wend = width;
%         return;
%     end
%     wstart = max(idx(1) - 4, 1);
%     wend = min(idx(end) + 4, width);
% 
%     
%     idx = [];
%     h_avg = mean(binary, 2);
%     idx = find(h_avg < 0.98);
%     if isempty(idx)
%         res = binary;
%         position.hstart = 0;
%         position.hend = height;
%         position.wstart = 0;
%         position.wend = width;
%         return;
%     end
%     hstart = max(idx(1) - 4, 1);
%     hend = min(idx(end) + 4, height);

    
%     idx = [];
    tmp = ones(length((1+step):(height-step)), 1);
    for i = (1+step):(width-step)
        tmp = tmp & binary((1+step):(height-step), i);
    end
    idx = find(tmp==0);
    if isempty(idx)
        res = binary;
        return;
    end
    hstart = max(step + idx(1) - 6, 1);
    hend = min(step + idx(end) + 12, height);
    
    idx = [];
    tmp = ones(1, length((1+step):(width-step)));
    for i = (1+step):(height-step)
        tmp = tmp & binary(i, (1+step):(width-step));
    end
    idx = find(tmp==0);
    if isempty(idx)
        res = binary;
        return;
    end
    wstart = max(step + idx(1) - 4, 1);
    wend = min(step + idx(end) + 4, width);
    
    res = [];
    res = img(hstart:hend, wstart:wend);
%     position.hstart = hstart;
%     position.hend = hend;
%     position.wstart = wstart;
%     position.wend = wend;
%     
end

