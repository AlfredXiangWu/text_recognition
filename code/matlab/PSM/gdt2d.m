function [varargout] = gdt2d(varargin)

% gdt2d:  2D generalized distance transform
%
% dt = gdt2d(M)

fprintf('Mex file not found -- attempting to compile...\n');
try
    mex gdt2d.cpp
    fprintf('Compilation complete.  Continuing.\n');
    [varargout{1:nargout}] = gdt2d(varargin{:});
catch
    error('Unable to run mex file.');
end
