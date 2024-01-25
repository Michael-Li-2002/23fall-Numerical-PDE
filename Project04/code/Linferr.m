function Linf_err = Linferr(pde,uh,xmesh,t)
%LINFERR L^infty error
N = size(uh,2); h = xmesh(2)-xmesh(1);
u = pde.exactu([t*ones(N,1),xmesh']);
Linf_err = max(abs(u-(uh)'));

