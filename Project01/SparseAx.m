function [V] = SparseAx(A,U,rN,tN)
%SPARSEAX AU=V
% Only for hw2. For A sparse and U vector, calculate AU=V
V = zeros(rN-1,tN);
for j = 1:tN
    V(1,j) = A(j,tN+j) * U(2,j);
    V(1,j) = V(1,j) + sum(A(j,1:tN).* U(1,1:tN));
end
for i = 2:rN-2
    j = 1; index = (i-1)*tN + j;
    V(i,j) = A(index,index)*U(i,j) + A(index,index+1)*U(i,j+1) + A(index,index+tN-1)*U(i,j+tN-1) + A(index,index+tN)*U(i+1,j) + A(index,index-tN)*U(i-1,j);
    for j = 2:tN-1
        index = (i-1)*tN + j;
        V(i,j) = A(index,index)*U(i,j) + A(index,index+1)*U(i,j+1) + A(index,index-1)*U(i,j-1) + A(index,index+tN)*U(i+1,j) + A(index,index-tN)*U(i-1,j);
    end
    j = tN; index = (i-1)*tN + j;
    V(i,j) = A(index,index)*U(i,j) + A(index,index+1-tN)*U(i,j+1-tN) + A(index,index-1)*U(i,j-1) + A(index,index+tN)*U(i+1,j) + A(index,index-tN)*U(i-1,j);
end
i = rN-1;
j = 1; index = (i-1)*tN + j;
V(i,1) = A(index,index)*U(i,j) + A(index,index+1)*U(i,j+1) + A(index,index+tN-1)*U(i,j+tN-1) + A(index,index-tN)*U(i-1,j);
for j = 2:tN-1
    index = (i-1)*tN + j;
    V(i,j) = A(index,index)*U(i,j) + A(index,index+1)*U(i,j+1) + A(index,index-1)*U(i,j-1) + A(index,index-tN)*U(i-1,j);
end
j = tN; index = (i-1)*tN + j;
V(i,j) = A(index,index)*U(i,j) + A(index,index+1-tN)*U(i,j+1-tN) + A(index,index-1)*U(i,j-1) + A(index,index-tN)*U(i-1,j);
end
