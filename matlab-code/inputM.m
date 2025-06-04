function data_x1 = inputM(age)
% Input the single cell data.
% Change the 0 row in data_x to avoid some errors in iterations.
% OUTPUT:
%   data_x1: the matrix will be used in the matrix factorization.

filename = sprintf('%s%d%s','..\data\origin\hv',age,'.csv');
file = importdata(filename);
data_x = file.data;
data_x1 = data_x';
if age == 19
    data_x1(:,all(data_x1==0,1)) = []; % delete all 0 column.
else
    gene = importdata('..\data\origin\age19gene0.csv');
    data_x1(:,gene) = [];
end
end
