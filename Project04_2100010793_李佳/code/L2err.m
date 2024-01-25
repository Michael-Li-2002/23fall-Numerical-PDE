function [L2_err] = L2err(pde,uh,xmesh,t)
%L2ERR L^2 error
N = size(uh,2); h = xmesh(2)-xmesh(1);
u = pde.exactu([t*ones(N,1),xmesh']);
L2_err = sqrt(h*sum((u-(uh)').^2));

