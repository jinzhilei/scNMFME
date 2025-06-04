function [a19, b19, t19] = pselect(data19,rank,para)
% Select alpha, beta and theta.
% INPUTS:
%   data19: transposed gene expression matrix at 19-year-olds (rows are cells and columns are genes)
%   rank: the selected cluster number
%   para: range of parameters
% OUTPUTS:
%   a19: $\alpha$ will be used in the factorization
%   b19: $\beta$ will be used in the factorization
%   t19: $\theta$ will be used in the factorization

alpha = para.alpha;
beta = para.beta;
theta = para.theta;
max_it = para.max_it;
lambda0 = 0; % for this selection using 19-year-old matrix, lambda is a useless parameter.
na = length(alpha);
nb = length(beta);
nt = length(theta);
times = 10;
abtef = zeros(1000,5); % 1000 should be a number that is big enough to contain all combinations.
nr = 0; % nr is the number of combination.
for i=1:nt
    for j = 1:nb
        for k = 1:na
            setdemorandstream(12345);
            ec = 0; % ec is the cumulative value of error.
            fc = 0;
            for u = 1:times
                whef = nmf(data19,rank,alpha(k),beta(j),theta(i),lambda0);
                ec = ec + mean(whef{1,3}((max_it-1):max_it));
                fc = fc + whef{1,4};
            end
            nr = nr + 1;
            abtef(nr,1:3) = [alpha(k),beta(j),theta(i)];
            abtef(nr,4) = ec/times;
            abtef(nr,5) = fc/times;
        end
    end
end
abtef(all(abtef==0,2),:) = [];
dlmwrite('..\data\output\abt_pselect.csv',abtef,',');

min_col4 = min(abtef(:,4));
row_col4 = find(abtef(:,4) == min_col4);
min_col5 = min(abtef(row_col4,5));
result_row = row_col4(find(abtef(row_col4,5) == min_col5));
[a19, b19, t19, ~, ~] = abtef(result_row, :);
end