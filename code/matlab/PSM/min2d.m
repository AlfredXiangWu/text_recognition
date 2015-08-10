%  MIN2D  Get the minimum element and indices from a 2d matrix
%
%  function [m,i,j] = min2d(M)

function [m,i,j] = min2d(M)

[m1,i1] = min(M,[],1);
[m,j] = min(m1,[],2);
i = i1(j);
