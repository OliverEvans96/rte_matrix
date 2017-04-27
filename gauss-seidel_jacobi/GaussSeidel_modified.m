% Sources:
% https://www.mathworks.com/matlabcentral/fileexchange/32051-gauss-seidel-method?focused=5193371&tab=function
tic
format compact

n = length(b);

% X is initial guess
X = zeros(n,1);
Error_eval = ones(n,1);
delta = 0.001;
max1=100;

iteration = 0;
while max(Error_eval) > delta
    iteration = iteration + 1
    Z = X;  % save current values to calculate error later
    for i = 1:n
        j = 1:n; % define an array of the coefficients' elements
        j(i) = [];  % eliminate the unknow's coefficient from the remaining coefficients
        Xtemp = X;  % copy the unknows to a new variable
        Xtemp(i) = [];  % eliminate the unknown under question from the set of values
        X(i) = (b(i) - sum(A(i,j) * Xtemp)) / A(i,i);
    end
    Xsolution(:,iteration) = X;
    Error_eval = sqrt((X - Z).^2);
    max(Error_eval);
    if iteration>=max1
        %break if at maximum number of iterations
     break
    end
end
toc