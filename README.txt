# scNMFME

----------
1. Overview:
This is a MATLAB and R package of scNMFME ("single-cell Non-negative Matrix Factorization Module Extraction"). 
scNMFME is a non-negative matrix factorization framework for quantifying immunosenescence from single-cell RNA-seq data.
----------

----------
2. Systems Requirements:
Both MATLAB and R are required for running scNMFME.
This Package has been tested using MATLAB 2023a and R 4.1.3 (with RStudio 2023.09.1) on Windows (64-bit).
----------

----------
3. Usage:

3.1 Main scripts in File "matlab-code":
(1) loopfunc.m -- main function including complete NMF process
(2) inputM.m -- input the single cell data and delete the 0 row in the matrix to avoid errors
(3) gap_cluster.m -- compute the optimal cluster number in the single cell matrix
(4) errorCalc.m -- calculate the error
(5) paraSet.m -- set the range of parameters
(6) pselect.m -- select the optimal values of alpha, beta, and theta
(7) lselect.m -- select the optimal value of lambda
(8) nmf.m -- the NMF function
(9) updateW.m -- update the matrix "W" (the module expression matrix)
(10) updateH.m -- update the matrix "H" (the module-gene matrix) for age 19
(11) updateH2.m -- update the matrix "H" (the module-gene matrix) for age 20-97

3.2 Main scripts in File "R-code":
Four scripts that need to be run in R. 
(1) data-preprocessing.R; 
(2) module-gene-matrix.R; 
(3) multi-age-distribution.R; 
(4) signature-scores.R
----------