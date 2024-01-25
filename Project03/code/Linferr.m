function Linf_err = Linferr(pde,uh,t)
%LINFERR L^infty error
N = size(uh,2)-1; h = 1/N;
u = pde.exactu([t*ones(N+1,1),h*(0:N)']);
Linf_err = max(abs(u-(uh)'));

