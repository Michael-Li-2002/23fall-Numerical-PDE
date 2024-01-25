%% Project 03 ��� 2100010793
% special case for \theta = 1/2, 1/2-1/(12\mu) and data_2 or data_3
% when using nonconstant time step, they will have better performance
clear all; close all;
%% README
% t_0 Ϊĥ��ʱ��, t_0ǰ��С����, t_0���ô󲽳�
% ÿ�����ж�����N����, ��Ҫ�� Parameter �и���:
% һ����ֹʱ��t, ÿ��N��Ӧ��\tau(��Ҫ��t_0ǰ��ֱ��), \theta(��t_0ǰ�� ����ֱ��ע�ͷ�ע��),
% ���к����N, ������, L^\infty, L^2 ���ı��, 
% ������ L^2 ���ֱ��� N ����������˫����ͼ
%% Parameters
option.Nlist = [8;16;32;64;128];                 % space number
option.hlist = 1./option.Nlist;                  % space step
option.t = 10;                                   % final time
option.t0 = 0.005;                               % time for smoothing initial data
option.taulist1 = 5/6 *option.hlist.^2   ;       % time step(before t_0)
option.taulist2 = 1 *option.hlist.^1   ;         % time step(after  t_0)
option.mulist1 = option.taulist1./option.hlist.^2; % mesh ratio
option.mulist2 = option.taulist2./option.hlist.^2; % mesh ratio
% option.theta1 = 0.5 *ones(size(option.Nlist,1),1);  % theta ( when theta selected const )
% option.theta2 = 0.5 *ones(size(option.Nlist,1),1);  % theta ( when theta selected const )
option.theta1 = 0.5 - 1./(12*option.mulist1);    % when theta = 1/2 - 1/(12\mu)
option.theta2 = 0.5 - 1./(12*option.mulist2);    % when theta = 1/2 - 1/(12\mu)
option.fds = @theta_HeatEq;                      % finite difference method: theta scheme
pde = data_2;                                    % choose initial data
                                                 % 1:smooth  2:continuous  3:piecewise continuous

%% 
num = size(option.Nlist,1);
Linf_err = zeros(num,1);
L2_err = zeros(num,1);
totalcal = zeros(num,1);

for i = 1:num
    N = option.Nlist(i); h = 1/N; 
    tau1 = option.taulist1(i); mu1 = tau1/h^2;
    tau2 = option.taulist2(i); mu2 = tau2/h^2;
    theta1 = option.theta1(i); theta2 = option.theta2(i);
    xmesh = h * (0:N)';
    M1 = ceil(option.t0/tau1);
    M2 = ceil((option.t-M1*tau1)/tau2);
    M = M1+M2; real_t = M1*tau1 + M2*tau2; % real final time
    u = zeros(M,N+1);

    %% Solve Heat Equation & Record calculations
    u0 = (pde.initdata(xmesh))';
    [u(1:M1,:),totalcal1] = option.fds(theta1,u0,mu1,M1);
    [u(M1+1:M,:),totalcal2] = option.fds(theta2,u(M1,:),mu2,M2);
    totalcal(i) = totalcal1 + totalcal2;
    Linf_err(i) = Linferr(pde,u(M,:),real_t);
    L2_err(i) = L2err(pde,u(M,:),real_t);
end

%% Display data & figure
fprintf('%4s | %11s | %8s | % 6s\n','N','calculation','Linf err','L2 err')
fprintf('%s\n',repmat('-',1,40))
fprintf('%4d | %11d | %8.2e | %6.2e\n',[(option.Nlist)'; (totalcal)'; (Linf_err)'; (L2_err)'])
fprintf('\n')

subplot(1,2,1);
loglog(option.Nlist,L2_err)
subplot(1,2,2);
loglog(totalcal,L2_err)