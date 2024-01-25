%% NMPDE Project02 上机作业
% 在Data Setteings中可改变参数\epsilon的大小
% 在Setting Nonuniform Mesh中可设置边界点bd, 及两侧的剖分数M,N
% 直接单击运行可得到数值解与真解的图像, 并输出无穷范数误差, 可复现上机报告的结果
%% Data Settings 
eps = 0.001;
u = @(x) ( 1 - exp(-x/eps) ) / ( 1 - exp(-1/eps) );
f = @(x) 1/eps^2 * exp(-x/eps) / ( 1 - exp(-1/eps) );
fDeriv = @(x) -1/eps^3 * exp(-x/eps) / ( 1 - exp(-1/eps) );
%% Setting Nonuniform Mesh
bd = 2 * eps * log(1/eps);
N = 16; M = 20;
% bd = 0.5;       % If uniform mesh is needed, uncomment these lines. 
% N = 40; M = N;  % Change the value of N, you can get uniform mesh of n=2N nodes.
n = N+M;
hN = (1 - bd)/N;
hM = bd/M;

%% Setting Linear System
A = zeros(n-1);
F = zeros(n-1,1); U = zeros(n-1,1);
% the coefficient matrix
diagline = ones(M-1,1);
A(1:M-1,1:M-1) = full(spdiags([ -1/hM^2 * diagline, 2/hM^2 * diagline, -1/hM^2 * diagline ], [-1,0,1], M-1,M-1));
A(M-1,M) = -1/hM^2; 
diagline = ones(N-1,1);
A(M+1:M+N-1,M+1:M+N-1) = full(spdiags([ -1/(hM*hN) * diagline, 2/(hM*hN) * diagline, -1/(hM*hN) * diagline ], [-1,0,1], N-1,N-1));
A(M+1,M) = -1/(hM*hN);
A(M,M-1) = -1/hM^2; A(M,M) = (hM+hN)/(hM^2 * hN); A(M,M+1) = -1/(hM*hN);

% Right hand term
for i = 1:M-1
    F(i) = f(i*hM);
end
F(M) = f(M*hM) - fDeriv(M*hM) * (hN-hM)/3;
for j = 1:N-1
    F(M+j) = f(bd + j*hN) * hN/hM;
end
F(n-1) = F(n-1) + 1/(hM*hN);

%% Solving Linear System Using CG
% initial value: 0
r0 = F; p = r0; 
alpha = 0; beta = 0; k = 0;
while max(max(abs(r0))) > 0.0001
    V = A*p;
    alpha = sum(sum(r0.*r0))/sum(sum(p.*V));
    U = U + alpha * p;
    r1 = r0 - alpha * V;
    beta = sum(sum(r1.*r1))/sum(sum(r0.*r0));
    r0 = r1;
    p = r0 + beta * p;
    k = k + 1;
    disp(['after ',num2str(k),' steps, the infty norm residue is ',num2str(max(max(abs(r0))))])
end

%% Plot and Print Numerical Results 
ax = [hM*(0:M),bd + hN*(1:N)]; % x-axis
UData = [0;U;1];               % Numerical solution
UReal = u(ax);                 % Real solution
plot(ax,UReal);
hold on;
plot(ax,UData);
% print the error in infinity norm
disp(['L^infty error : ',num2str(max(abs(UData-UReal')))]);