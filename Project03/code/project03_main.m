%% Project 03 李佳 2100010793
% Only for simple case: constant time step; 
% when theta\neq 0,1, \tau=O(h) in data_2 or data_3,
% oscillation may occur. Project_specialcase deal with them.
%% README
% 每次运行对所有N计算, 需要在 Parameter 中给定:
% 一个终止时间t, 每个N对应的\tau, 一个\theta,
% 运行后输出N, 运算量, L^\infty, L^2 误差的表格, 
% 并作出 L^2 误差分别与 N 、运算量的双对数图
%% Parameters
clear all; close all;
option.Nlist = [8;16;32;64;128];                 % space number
option.hlist = 1./option.Nlist;                  % space step
option.t = 1;                                    % final time
option.taulist = 1 *option.hlist.^1   ;          % time step(change the coefficient and power) 
option.mulist = option.taulist./option.hlist.^2; % mesh ratio
option.theta = 1/2*ones(size(option.Nlist,1),1); % theta ( when theta selected const )
% option.theta = 0.5 - 1./(12*option.mulist);    % when theta = 1/2 - 1/(12\mu)
option.fds = @theta_HeatEq;                      % finite difference method: theta scheme
pde = data_2;                                    % choose initial data
                                                 % 1:smooth  2:continuous  3:piecewise continuous

%% 
num = size(option.Nlist,1);  
Linf_err = zeros(num,1);         % L^\infty error 
L2_err = zeros(num,1);           % L^2 error 
totalcal = zeros(num,1);         % total calculation

for i = 1:num
    N = option.Nlist(i); h = 1/N; 
    tau = option.taulist(i); 
    mu = tau/h^2;
    theta = option.theta(i);
    xmesh = h * (0:N)';      % spatial uniform mesh                    
    M = ceil(option.t/tau);  
    real_t = M * tau;        % real final time
    u = zeros(M,N+1);

    %% Solve Heat Equation & Record calculations
    u0 = (pde.initdata(xmesh))';                  % initial data
    [u,totalcal(i)] = option.fds(theta,u0,mu,M);  
    Linf_err(i) = Linferr(pde,u(M,:),real_t);
    L2_err(i) = L2err(pde,u(M,:),real_t);
end

%% Display data & figure
fprintf('%4s | %11s | %8s | % 6s\n','N','calculation','Linf err','L2 err')
fprintf('%s\n',repmat('-',1,40))
fprintf('%4d | %11d | %8.2e | %6.2e\n',[(option.Nlist)'; (totalcal)'; (Linf_err)'; (L2_err)'])
fprintf('\n')

L2data = L2_err';
caldata = totalcal';
subplot(1,2,1);
loglog(option.Nlist,L2_err)
subplot(1,2,2);
loglog(totalcal,L2_err)