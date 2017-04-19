% Sources:
% http://www.mathworks.com/matlabcentral/fileexchange/2182-numerical-methods-using-matlab--3e?focused=5042682&tab=function

tic;
X=0;
N = length(b);
% P is an N x 1 matrix; the initial guess
P = zeros([N,1]);
% delta is the tolerance for P
delta = 0.001;
% max 1 is the maximum number of iterations
max1 = 10;

%for the set number of maximum number of iterations
for k=1:max1
   for j=1:N
      % Develop each row of x^(i+1) as ((b-(U+L))*x^i)/D
      X(j)=(b(j)-A(j,[1:j-1,j+1:N])*P([1:j-1,j+1:N]))/A(j,j);
   end
   %calculate error and relative error
   err=abs(norm(X'-P))
   relerr=err/(norm(X)+eps);
   %set new x^i
   P=X';
    %break from loop if error is within tolerance
      if (err<delta)|(relerr<delta)
     break
   end
end

X=X';
toc;