%% Project 03 ¿Óº— 2100010793
% special case for \theta = 1/2, 1/2-1/(12\mu) and data_3
% when using nonconstant time step, they will have better performance
clear all; close all;
%% Parameters
option.Nlist = [8;16;32;64;128];                 % space number
option.hlist = 1./option.Nlist;                  % space step
option.t0 = 0.005;                                 % time for smoothing initial data
option.fds = @theta_HeatEq;                      % finite difference method: theta scheme
pde = data_3;                                    % choose initial data
                                                 % 1:smooth  2:continuous  3:piecewise continuous

option.taulist1 = [option.hlist.^2,5/6 * option.hlist.^2]   ;
option.taulist2 = [1 *option.hlist.^1,1 *option.hlist.^1]   ;            % time step ()
option.mulist1 = option.taulist1./option.hlist.^2;
option.mulist2 = option.taulist2./option.hlist.^2; % mesh ratio
% option.theta1 = 0.5 *ones(size(option.Nlist,1),1);  % theta ( when theta selected const )
% option.theta2 = 0.5 *ones(size(option.Nlist,1),1);  % theta ( when theta selected const )
option.thetalist1 = [0.5*ones(size(option.Nlist,1),1),0.5 - 1./(12*option.mulist1)];     % when theta = 1/2 - 1/(12\mu)
option.thetalist2 = [0.5*ones(size(option.Nlist,1),1),0.5 - 1./(12*option.mulist2)];     % when theta = 1/2 - 1/(12\mu)
option.tlist = [0.1,1,10];

%%
L2data = zeros(6,5);
caldata = zeros(6,5);
for ind = 1:3
    t = option.tlist(ind);
    for k=1:2
        taulist1 = option.taulist1(:,k); taulist2 = option.taulist2(:,k);            % time step ()
        mulist1 = option.mulist1(:,k);  mulist2 = option.mulist2(:,k);
        thetalist1 = option.thetalist1(:,k);  % theta ( when theta selected const )
        thetalist2 = option.thetalist2(:,k);
        %% 
        num = size(option.Nlist,1);  
        Linf_err = zeros(num,1);         % L^\infty error 
        L2_err = zeros(num,1);
        totalcal = zeros(num,1);         % total calculation

        for i = 1:num
            N = option.Nlist(i); h = 1/N; 
            tau1 = taulist1(i); mu1 = tau1/h^2;
            tau2 = taulist2(i); mu2 = tau2/h^2;
            theta1 = thetalist1(i); theta2 =thetalist2(i);
            xmesh = h * (0:N)';
            M1 = ceil(option.t0/tau1);
            M2 = ceil((t-M1*tau1)/tau2);
            M = M1+M2; real_t = M1*tau1 + M2*tau2;
            u = zeros(M,N+1);

            %% Solve Heat Equation & Record calculations
            u0 = (pde.initdata(xmesh))';
            [u(1:M1,:),totalcal1] = option.fds(theta1,u0,mu1,M1);
            % u(M1,:) = (pde.exactu([M1*tau1*ones(N+1,1),xmesh]))';
            [u(M1+1:M,:),totalcal2] = option.fds(theta2,u(M1,:),mu2,M2);
            totalcal(i) = totalcal1 + totalcal2;
            Linf_err(i) = Linferr(pde,u(M,:),real_t);
            L2_err(i) = L2err(pde,u(M,:),real_t);
        end
        Nedata = [option.Nlist,L2_err];
        Cedata = [totalcal,L2_err];
        filename1 = ['Ne_3',num2str(ind),num2str(k+5),'.dat'];
        filename2 = ['Ce_3',num2str(ind),num2str(k+5),'.dat'];
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
        disp(['finish ',num2str((ind-1)*2+k),'/',num2str(6)])
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