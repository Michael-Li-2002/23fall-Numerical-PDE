function [uh,totalcal] = theta_HeatEq(theta,u0,mu,M)
%
N = size(u0,2)-1;
u = zeros(M+1,N+1);
u(1,:) = u0;
totalcal = 0;

% coeffcient matrix of U^{m+1}
A = (1+2*theta*mu)*eye(N-1) + sparse([(2:N-1)';(1:N-2)'],[(1:N-2)';(2:N-1)'],...
    -theta*mu*ones(2*(N-2),1),N-1,N-1);
if theta == 0
    for m = 1:M
        u(1+m,2:N) = mu * (u(m,1:N-1) + u(m,3:N+1)) + (1-2*mu) * u(m,2:N);
        totalcal = totalcal + (N-1)*5; % 3 time and 2 plus
    end
elseif theta == 1
    for m = 1:M
        [u(1+m,2:N),cal] = tridiagSolver(A,u(m,2:N));
        totalcal = totalcal + cal;
    end
else
    for m = 1:M
        mid = (1-theta)*mu * (u(m,1:N-1) + u(m,3:N+1)) + (1-2*(1-theta)*mu) * u(m,2:N);
        totalcal = totalcal + (N-1)*5; % 3 time and 2 plus
        [temp,cal] = tridiagSolver(A,(mid)');
        u(1+m,2:N) = temp';
        totalcal = totalcal + cal;
    end 
end
uh = u(2:M+1,:);
