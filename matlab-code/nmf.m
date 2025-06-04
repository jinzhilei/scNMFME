function whef = nmf(age,rank,alpha,beta,theta,lambda)
% This is the NMF function.
% The objective function: (||X-WH||)^2 + alpha*||RH|| + beta*(\sum cos(W_i,W_j) + theta*(||W||_2,1/2)^1/2.
% INPUTS:
%   rank: the selected cluster number
%   para: range of parameters
%   alpha, beta, theta, lambda: regularization parameters
% OUTPUT:
%   whef = {W, H, E, F}
%   W: the module expression matrix
%   H: the gene-module matrix
%   E: Error = (||X-WH||)^2
%   F: the objective function value

data_x = inputM(age);
para = paraSet();
max_it = para.max_it;
error1 = para.error;

nx = size(data_x);
H = 10 * rand(rank, nx(2));
H = H ./sum(H,2);
W = 10 * rand(nx(1), rank);
W = W ./ sum(W,1);
E = zeros(max_it, 1);
if age == 19
    for i = 1:max_it 
        E(i) = errorCalc(data_x,W,H);
        if E(i) < error1
            break
        else
            WQS = updateW(data_x,W,H,theta,beta);
            W = WQS{1,1};
            HR = updateH(data_x,W,H,alpha);
            H = HR{1,1};
        end
    end
    Q = WQS{1,2};
    S = WQS{1,3};
    R = HR{1,2};
    F = trace(data_x' * data_x) - 2*trace(data_x' * W * H) + trace(H' * W' * W * H) + ...
        + alpha * trace(H' * R' * R * H) + 4*theta*trace(W' * Q * W) + beta * trace(W' * W * S);
whef = {W, H, E, F};
else
    H19 = importdata('..\data\output\h19_p00100110.csv');
    % H19 = importdata('..\GSE157007\outdata\h19_sam1.csv');
    for i = 1:max_it
        E(i) = errorCalc(data_x,W,H);
        if E(i) < error1
            break
        else
            WQS = updateW(data_x,W,H,theta,beta);
            W = WQS{1,1};
            size(W)
            HR = updateH2(data_x,W,H,H19,alpha,lambda);
            H = HR{1,1};
        end
    end
    Q = WQS{1,2};
    S = WQS{1,3};
    R = HR{1,2};
    F = trace(data_x' * data_x) - 2*trace(data_x' * W * H) + trace(H' * W' * W * H) + ...
        + alpha * trace(H' * R' * R * H) + 4*theta*trace(W' * Q * W) + beta * trace(W' * W * S) + ...
        + lambda * trace((H'-H19')*(H-H19));
whef = {W, H, E, F};
end
end