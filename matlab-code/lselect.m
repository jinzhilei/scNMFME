function l20 = lselect(data20,rank,alpha,beta,theta,para)
% Select lambda.
% INPUTS:
%   data20: transposed gene expression matrix at 20-year-olds (rows are cells and columns are genes)
%   rank: the selected cluster number
%   alpha, beta, theta: the selected parameters before
%   para: range of parameters
% OUTPUTS:
%   l20: $\lambda$ will be used in the factorization

max_it = para.max_it;
lambda = para.lambda;
nl = length(lambda);
times = 10;
lef = zeros(nl,3); % 1000 should be a number that is big enough to contain all combinations.
nr = 0; % nr is the number of combination.
for i=1:nl
    setdemorandstream(12345);
    ec = 0; % ec is the cumulative value of error.
    fc = 0;
    for u = 1:times
        whef = nmf(data20,rank,alpha,beta,theta,lambda(i));
        ec = ec + mean(whef{1,3}((max_it-1):max_it));
        fc = fc + whef{1,4};
    end
    nr = nr + 1;
    lef(nr,1) = lambda(i);
    lef(nr,2) = ec/times;
    lef(nr,3) = fc/times;
end
lef(all(lef==0,2),:) = [];
dlmwrite('..\data\output\l_select.csv',lef,',');

min_col2 = min(lef(:,2));
row_col2 = find(lef(:,2) == min_col2);
min_col3 = min(lef(row_col2,3));
result_row = row_col2(find(lef(row_col2,3) == min_col3));
[l20, ~, ~] = lef(result_row, :);
end