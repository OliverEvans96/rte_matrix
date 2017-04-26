function [egs,ejacobi] = Convergence_Check(A_mat)

% Sources:
% http://www.codewithc.com/gauss-seidel-method-matlab-program/


% Separation of matrix A into lower triangular and upper triangular matrices
% A = D + L + U
D = diag(diag(A_mat));
L = tril(A_mat)- D;
U = triu(A_mat)- D;

% check for convergence condition for Gauss-Seidel method
egs= max(abs(eigs(-inv(D+L)*(U))));
%if abs(egs) >= 1
%    disp ('Since the modulus of the largest Eigen value of iterative matrix is not less than 1') 
%    disp ('this process is not convergent for Gauss-Seidel.')
%end
%
% check for convergence condition for Jacobi method
ejacobi= max(abs(eigs(-inv(D)*(U+L))));
%if abs(ejacobi) >= 1
%    disp ('Since the modulus of the largest Eigen value of iterative matrix is not less than 1') 
%    disp ('this process is not convergent for Jacobi.')
%end

