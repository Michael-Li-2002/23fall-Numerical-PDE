%% Project 03 李佳 2100010793
% Only for simple case: constant time step; 
% when theta\neq 0,1, \tau=O(h) in data_3,
% oscillation may occur. Project_specialcase deal with them.
%% README
% 每次运行对所有N计算, 需要在 Parameter 中给定一个终止时间t, 
% 一个\theta,每个N对应的\tau. 运行后输出N, 运算量, L^\infty, L^2 误差的表格, 
% 并作出 L^\infty 误差分别与 N 、运算量的双对数图
%% Parameters
clear all; close all;
option.Nlist = [8;16;32;64;128];                 % space number
option.hlist = 1./option.Nlist;                  % space step
option.fds = @theta_HeatEq;                      % finite difference method: theta scheme
pde = data_3;                                    % choose initial data
                                                 % 1:smooth  2:continuous  3:piecewise continuous
option.taulist = [1/2 *option.hlist.^2, option.hlist.^1, 3*option.hlist.^2,option.hlist.^1, option.hlist.^2];
option.mulist = option.taulist./option.hlist.^2; % mesh ratio
option.thetalist = [0 *ones(5,1), 1 *ones(5,1),1 *ones(5,1),1/2 *ones(5,1),0.5 - 1./(12*option.mulist(:,5))];
option.tlist = [0.1,1,10];

L2data = zeros(15,5);
caldata = zeros(15,5);
for ind = 1:3
    t = option.tlist(ind);
    for k=1:5
        taulist = option.taulist(:,k);            % time step ()
        mulist = option.mulist(:,k);
        thetalist = option.thetalist(:,k);  % theta ( when theta selected const )

        %% 
        num = size(option.Nlist,1);  
        Linf_err = zeros(num,1);         % L^\infty error 
        L2_err = zeros(num,1);
        totalcal = zeros(num,1);         % total calculation

        for i = 1:num
            N = option.Nlist(i); h = 1/N; 
            tau = taulist(i); mu = tau/h^2;
            theta = thetalist(i);
            xmesh = h * (0:N)';
            M = ceil(t/tau); real_t = M * tau;     % real final time
            u = zeros(M,N+1);

            %% Solve Heat Equation & Record calculations
            u0 = (pde.initdata(xmesh))';                  % initial data
            [u,totalcal(i)] = option.fds(theta,u0,mu,M);  
            Linf_err(i) = Linferr(pde,u(M,:),real_t);
            L2_err(i) = L2err(pde,u(M,:),real_t);
        end
        Nedata = [option.Nlist,L2_err];
        Cedata = [totalcal,L2_err];
        filename1 = ['Ne_3',num2str(ind),num2str(k),'.dat'];
        filename2 = ['Ce_3',num2str(ind),num2str(k),'.dat'];
        fid = fopen(filename1,'a');
        for j=1:num
            fprintf(fid,'%d %e\n',Nedata(j,1),Nedata(j,2));
        end
        fclose(fid);
        fid = fopen(filename2,'a');
        for j=1:num
            fprintf(fid,'%d %e\n',Cedata(j,1),Cedata(j,2));
        end
        fclose(fid);
        disp(['finish ',num2str((ind-1)*5+k),'/',num2str(15)])
    end
end
%% Display data & figure

% fprintf('%4s | %11s | %8s | % 6s\n','N','calculation','Linf err','L2 err')
% fprintf('%s\n',repmat('-',1,40))
% fprintf('%4d | %11d | %8.2e | %6.2e\n',[(option.Nlist)'; (totalcal)'; (Linf_err)'; (L2_err)'])
% fprintf('\n')
% 
% L2data = L2_err';
% caldata = totalcal';
% subplot(1,2,1);
% loglog(option.Nlist,L2_err)
% subplot(1,2,2);
% loglog(totalcal,L2_err)