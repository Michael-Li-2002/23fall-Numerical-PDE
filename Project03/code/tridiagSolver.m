function [b,totalcal] = tridiagSolver(A,b)
% tridiagonal matrix linear equation solver

N = size(A,1);
totalcal = 0;
for i = 2:N
    b(i) = b(i) - b(i-1)*A(i,i-1)/A(i-1,i-1);
    totalcal = totalcal + 3;
    A(i,i) = A(i,i) - A(i-1,i)*(A(i,i-1)/A(i-1,i-1));
    totalcal  = totalcal + 3;
end
b(N) = b(N)/A(N,N);
totalcal = totalcal + 1;
for i = N-1:-1:1
    b(i) = b(i) - b(i+1)*A(i,i+1);
    totalcal = totalcal + 2;
    b(i) = b(i)/A(i,i);
    totalcal = totalcal + 1;
end
