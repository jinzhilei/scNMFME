function e = errorCalc(X,W,H)
% Calculate the error. 
% e = ||X-W*H||_F.
e = norm((X-W*H),'fro');
end