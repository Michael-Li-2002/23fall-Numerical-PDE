function [u,exactt] = Upwind_conservative(pde,xmesh,h,tau,t)
%UPWIND_CONSERVATIVE conservative upwind scheme

tstep = ceil(t/tau);
exactt = tau * tstep;
N = size(xmesh,2);
u = pde.initu(xmesh);
ind = 2:N-1;
a = zeros(1,N);
for i = 1: tstep
    for j = 1:N-1
        if u(j) ~= u(j+1)
            a(j) = ( pde.f(u(j+1))-pde.f(u(j)) )/ (u(j+1)-u(j));
        else
            a(j) = 0;
        end
    end
    u(ind) = u(ind) - tau/(2*h)*((1+sign(a(ind))).*pde.f(u(ind)) + (1-sign(a(ind))).*pde.f(u(ind+1))...
        -(1+sign(a(ind-1))).*pde.f(u(ind-1)) - (1-sign(a(ind-1))).*pde.f(u(ind))) ;
    u(1) = pde.lbd([tau*i]);
    u(N) = pde.rbd([tau*i]);
end

