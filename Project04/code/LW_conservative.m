function [u,exactt] = LW_conservative(pde,xmesh,h,tau,t)
%LW_CONSERVATIVE Conservative Lax-Wendroof scheme
tstep = ceil(t/tau);
exactt = tau * tstep;
N = size(xmesh,2);
u = pde.initu(xmesh);
ind = 2:N-1;
a = zeros(1,N);
for i = 1: tstep
    a(1:N-1) = pde.Df((u(1:N-1)+u(2:N))/2);
    u(ind) = u(ind) - tau/(2*h)*(pde.f(u(ind+1)) - pde.f(u(ind-1)) ) ...
        + tau^2/(2*h^2)*(a(ind).*(pde.f(u(ind+1)) - pde.f(u(ind)))  ...
        - a(ind-1).*(pde.f(u(ind)) - pde.f(u(ind-1))));
    u(1) = pde.lbd([tau*i]);
    u(N) = pde.rbd([tau*i]);
end

