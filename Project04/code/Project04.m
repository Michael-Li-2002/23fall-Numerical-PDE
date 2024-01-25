%% Project 04 2100010793 李佳
% 在设置中更改: 
%     option.h( x 方向步长 ), option.t( 终止时间 t ), option.tau( 时间步长)
%     option.FDS( 有限差分格式: Upwind, Upwind_conservative, LW_conservative )
% 可输出对应时刻的数值解图像及其与弱解的 L^\infty, L^2 误差(但由于解间断, L^\infty 误差的意义不大)
% 注释部分用于批量生成上机报告中的数据
%% Settings
option.xrange = [-1,1];
option.h = 10^(-3);                                % x direction step size
option.t = 1.0;                                    % final time
option.tau = 1 * option.h;                         % time step size
% option.hlist = [10^(-1),10^(-2),10^(-3),10^(-4)];
% option.tlist = 0.2* (1:5);
% option.taulist = option.hlist;
option.pde = Burgers();                            
option.FDS = @Upwind_conservative;                 % the finite difference schemes:
                                                   % Upwind, Upwind_conservative, LW_conservative
%% Compute numerical solution
N = ceil((option.xrange(2)-option.xrange(1))/option.h+0.1);
xmesh = option.xrange(1) + option.h * (0:N-1);
[uh,exactt] = option.FDS(option.pde,xmesh,option.h,option.tau,option.t);
Linf_err_disp = Linferr(option.pde,uh,xmesh,exactt);
L2_err_disp = L2err(option.pde,uh,xmesh,exactt);
disp(['L^\infty error: ',num2str(Linf_err_disp)])
disp(['L^2 error: ',num2str(L2_err_disp)])
plot(xmesh,uh);
%% Generate data
% Linf_err = zeros(size(option.hlist,2),size(option.tlist,2));
% L2_err = Linf_err;
% 
% for hind = 1:size(option.hlist,2)
%     h = option.hlist(hind); tau = option.taulist(hind);
%     N = ceil((option.xrange(2)-option.xrange(1))/h+0.1);
%     xmesh = option.xrange(1) + h * (0:N-1);
%     for tind = 1:size(option.tlist,2)
%         t = option.tlist(tind);
%         [uh,exactt] = option.FDS(option.pde,xmesh,h,tau,t);
%         Linf_err(hind,tind) = Linferr(option.pde,uh,xmesh,exactt);
%         L2_err(hind,tind) = L2err(option.pde,uh,xmesh,exactt);
%         filename = ['LWc',num2str(hind),num2str(tind),'.dat'];
%         fid = fopen(filename,'a');
%         fprintf(fid,'%e %e\n',[xmesh;uh]);
%         fclose(fid);
%     end
% end