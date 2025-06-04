function HR = updateH2(x1,W,H,H19,alpha,lambda)
% Update H for age20-97.
% INPUTS:
%   x1: the gene expression matrix
%   W: the module expression matrix after the previous iteration
%   H: the gene-module matrix after the previous iteration
%   H19: the result of gene expression matrix factorization at 19 years old 
%   alpha, lambda: regularization parameters
% OUTPUT:
%   HR = {H1,R} 

nh = size(H);
R = rand(nh(1), nh(1));
H1 = H .* ((W'*x1 + lambda*H19) ./ (W'*W*H + alpha*(R')*R*H + lambda*H + 1e-10));
HR = {H1,R};
end