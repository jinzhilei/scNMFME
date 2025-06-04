function para = paraSet()
% Set the parameters.
%   error1: One of the termination conditions for the matrix factorization iterative process
%   max_it: Maximum number of iterations for the matrix factorization iterative process
%   epsilon: the small constant is used to avoid overflow during computation
%   alpha, beta, theta, lambda are the range of regularization parameters
% The selection functions of regularization parameters are 'pselect' and 'lselect'.
para.error1 = 10;
para.max_it = 50;
para.epsilon = 1e-4;

para.alpha = [0.01,0.1,1,10];
para.beta = [0.01,0.1,1,10];
para.theta = [0.01,0.1,1,10];
para.lambda = [0.01,0.05,0.1,0.5,1,5,10,50,100];

end