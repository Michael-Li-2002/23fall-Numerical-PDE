function [L2_err] = L2err(pde,uh,t)
%L2ERR L^2 error
N = size(uh,2)-1; h = 1/N;
u = pde.exactu([t*ones(N+1,1),h*(0:N)']);
L2_err = sqrt(h*sum((u-(uh)').^2));

