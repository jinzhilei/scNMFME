function HR = updateH(x1,W,H,alpha)
% Update H for 19-year-old single cell data.
% INPUTS:
%   x1: the gene expression matrix
%   W: the module expression matrix after the previous iteration
%   H: the gene-module matrix after the previous iteration
%   alpha: regularization parameters
% OUTPUT:
%   HR = {H1,R} 

a = size(H);
R = rand(a(1), a(1));
H1 = H .* (W'*x1 ./ (W'*W*H + alpha*(R')*R*H + 1e-10));
HR = {H1,R};
end