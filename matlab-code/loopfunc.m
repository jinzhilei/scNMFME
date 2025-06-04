%% 19-97 NMF program.
% select k and parameters and implement matrix factorization 

%% Select k, alpha, beta and theta. Implement matrix factorization on 19-year-old gene expression matrix.
data19 = inputM(19);
rank = gap_cluster(data19);
para = paraSet();
[alpha,beta,theta] = pselect(data19,rank,para);
lambda0 = 0; % for 19-year-old matrix, lambda is a useless parameter.
whef19 = nmf(19,rank,alpha,beta,theta,lambda0);
W19 = whef19{1,1};
H19 = whef19{1,2};
dlmwrite('..\data\output\w19.csv', W19,',');
dlmwrite('..\data\output\h19.csv', H19,',');

%% Select lambda. Implement matrix factorization on 20-year-old gene expression matrix.
data20 = inputM(20);
lambda = lselect(data20,rank,alpha,beta,theta,para);
whef20 = nmf(20,rank,alpha,beta,theta,lambda);
W20 = whef20{1,1};
H20 = whef20{1,2};
dlmwrite('..\data\output\w20.csv', W20,',');
dlmwrite('..\data\output\h20.csv', H20,',');

%% Implement matrix factorization on 21-97 gene expression matrices.
for age = [21:91,93:97]
    data_x1 = inputM(age);
    whef = nmf(data_x1,rank,alpha,beta,theta,lambda);
    E = whef{1,3};
    if isnan(E(50))
        continue
    else
        W = whef{1,1};
        H = whef{1,2};
        wpath = sprintf('%s%d%s','..\data\output\w',age,'.csv');
        hpath = sprintf('%s%d%s','..\data\output\h',age,'.csv');
        dlmwrite(wpath, W,',');
        dlmwrite(hpath, H,',');
    end
end