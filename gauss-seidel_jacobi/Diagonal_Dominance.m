function [Aprecond,bprecond,Vprecond] = Diagonal_Dominance(A,b,nsv)
% Sources:
% Yuan, Jin Yun and Plamen Y. Yalamov.  A Method for Contructing Diagonally Dominant Preconditioners based on Jacobi Rotations.
% This code attempts to implement proposed algorithm to efficiently
% establish increased diagonal dominance

% Choose m (the number of Jacobi iterations to be applied), and ? > 0(the level of diagonal dominance);

% nsv = # of singular values to use for SVDS
tic
m=30000;
sigma=0.001;
A1=A;
B1=b;
N = length(B1);
Vprecond = eye(N);
for k=1:m
    %Find ai0j0 such that |ai0j0| = maxi6=j |aij |;
    D = diag(diag(A1));
    UL = A1-D;
	k
    [M,I]=max(abs(UL(:)))
    [i0, j0] = ind2sub(size(UL),I);
    % i1 = min(i0, j0); j1 = max(i0, j0);
    i1=min(i0,j0);
    j1=max(i0,j0);
    % Compute Qk =(ai1i1 ai1j1;aj1i1 aj1j1) = = UkSkVk^T (the SVD of matrix Qk);
    Qk=[A1(i1,i1),A1(i1,j1);A1(j1,i1),A1(j1,j1)];
    [U,S,V]=svds(Qk,nsv);
    %Apply the transformation Uk^T to rows i1 and j1 of matrix A
    Uk=eye(N);
    Uk(i1,i1)=U(1,1);
    Uk(i1,j1)=U(1,2);
    Uk(j1,i1)=U(2,1);
    Uk(j1,j1)=U(2,2);
    Ukt=transpose(Uk);
    A1=Ukt*A1;
    %Apply the transformation Vk to rows i1 and j1 of matrix A
    Vk=eye(N);
    Vk(i1,i1)=V(1,1);
    Vk(i1,j1)=V(1,2);
    Vk(j1,i1)=V(2,1);
    Vk(j1,j1)=V(2,2);
    Vprecond=Vprecond*Vk;
    %A1=Vk*A1;
	A1=A1*Vk;
    %Compute b = (U(k))T b;
    B1=Ukt*B1;
end
% Compute l (the number of rows in which |aii| ? sum |aij | ? ? );
l=0;
for i = 1:N
    j = 1:N;
    j(i) = [];
    B = abs(A1(i,j));
    check = abs(A1(i,i)) - sum(B); % Is the diagonal value greater than the remaining row values combined?
    if check >= sigma
        l=l+1;
    end
end
% Generate a preconditioner M: M(1 : l, 1 : l) = A(1 : l, 1 : l), M(i, i) = A(l, l), i = l + 1, . . . , n;
M1=A(l,l)*eye(N,N);
M1(1:l,1:l)=A1(1:l,1:l);
% Solve Ay = b with a preconditioner M;
Aprecond=M1*A1;
bprecond=M1*B1;
toc
