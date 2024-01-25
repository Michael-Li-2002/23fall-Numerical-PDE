function [u,exactt] = Upwind(pde,xmesh,h,tau,t)
%UPWIND Nonconservative Upwind scheme

tstep = ceil(t/tau);
exactt = tau * tstep;
N = size(xmesh,2);
u = pde.initu(xmesh);
ind = 2:N-1;
for i = 1: tstep
    uwind = ind - sign(pde.Df(u(ind)));
    u(ind) = (1-tau/h * abs(pde.Df(u(ind)))).*u(ind) + tau/h * abs(pde.Df(u(ind))).* u(uwind);
    u(1) = pde.lbd([tau*i]);
    u(N) = pde.rbd([tau*i]);
end

