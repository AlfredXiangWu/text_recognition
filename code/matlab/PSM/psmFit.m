function [varargout] = psmFit(varargin)

% psmFit:  Fits a part-structured model to a binary image.
%
% Uses bicubic translation for subpixel accuracy on translation.
%
% [dtsq,loc] = psmFit(model,img)
% [dtsq,loc] = psmFit(model,img,mdsq)
% Note:  mdsq = max_x^2+max_y^2 for equivalence with psmFit_gpu
%
% See also psmFit_gpu.

fprintf('Mex file not found -- attempting to compile...\n');
try
    mex psmFit.cpp
    fprintf('Compilation complete.  Continuing...\n');
    [varargout{1:nargout}] = psmFit(varargin{:});
catch
    error('Unable to run mex file.');
end
