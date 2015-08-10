function [varargout] = psmFit_gpu(varargin)

% psmFit_gpu:  Fits a part-structured model to a binary image, 
%   using parallel gpu code
%
% Uses bicubic translation for subpixel accuracy on translation.
%
% [dtsq,loc] = psmFit_gpu(model,img,[max_x max_y])
%
% See also psmFit.

fprintf('Mex file not found -- attempting to compile...\n');
try
    % These commands compile the CUDA source code on my system.
    % Yours may require something different depending upon how it is
    % configured.  I can't offer much support here, unfortunately.
    % Please check out nVidia's CUDA documentation for help.
    fname = 'psmFit_gpu';
    system(sprintf('nvcc -I"%s/extern/include" --cuda "%s.cu" --output-file "%s.cpp"',matlabroot,fname,fname));
    eval(sprintf('mex -L"C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v4.0/lib/x64" -lcufft -lcudart %s.cpp',fname));

    fprintf('Compilation complete.  Continuing...\n');
    [varargout{1:nargout}] = psmFit(varargin{:});
catch
    error('Unable to run mex file.');
end
