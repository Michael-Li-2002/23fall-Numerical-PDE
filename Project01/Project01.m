%% NMPDE Project 01 上机作业 李佳
% 在 Data Setteings 一节中可改变参数值\alpha,\beta,R_0的大小
% 在 Polar Coordinate Mesh 一节中可改变径向与角度上网格的均匀剖分数
% 请依次运行各节, 运行完成 Solving Linear System by CG 这一节后, 
% 依次运行后面的各节可依次复现上机报告的结果
% (y轴上Deflection的值为了方便和Load放在一起, 作图时作了比例调整. 该节中请手动调整比例)
%% Data Settings
alpha = 4; beta = 8; R0 = 0.6;
pressureXY = @(x,y) alpha * exp(-beta^2*(x^2+(y-R0)^2));
pressureRT = @(r,theta) pressureXY(r*cos(theta),r*sin(theta));

%% Polar Coordinate Mesh
rN = 160; rstep = 1/rN;    % r - radial direction
tN = 160; tstep = 2*pi/tN; % t - theta  direction
n = (rN-1)*tN;
U = zeros(rN-1,tN); % Deflection
P = zeros(rN-1,tN); % Load

%% Setting Linear System 
A = zeros(n); % Record the coefficient matrix
diagline = ones(tN,1); I = eye(tN);
T = full(spdiags([-(1/(tstep^2))*diagline,(+2+2/(tstep^2))*diagline,-(1/(tstep^2))*diagline],[-1, 0, 1],tN,tN));
T(1,tN) = -1/(tstep^2); T(tN,1) = -1/(tstep^2);
A(1:tN,1:tN) = T - 1/(2*rN) * ones(tN); A(1:tN,tN+1:2*tN) = -1.5*I;
for i = 2:rN-2
    T = full(spdiags([(-1/(i*tstep^2))*diagline,(2*i+2/(i*tstep^2))*diagline,(-1/(i*tstep^2))*diagline],[-1, 0, 1],tN,tN));
    T(1,tN) = -1/(i*tstep^2); T(tN,1) = -1/(i*tstep^2);
    A((i-1)*tN+1:i*tN,(i-1)*tN+1:i*tN) = T;
    A((i-1)*tN+1:i*tN,(i-2)*tN+1:(i-1)*tN) = -(2*i-1)/2*I;
    A((i-1)*tN+1:i*tN,(i)*tN+1:(i+1)*tN) = -(2*i+1)/2*I;
end
i = rN-1;
T = full(spdiags([-(1/(i*tstep^2))*diagline,(+2*i+2/(i*tstep^2))*diagline,-(1/(i*tstep^2))*diagline],[-1, 0, 1],tN,tN));
T(1,tN) = -1/(i*tstep^2); T(tN,1) = -1/(i*tstep^2);
A((i-1)*tN+1:i*tN,(i-1)*tN+1:i*tN) = T;
A((i-1)*tN+1:i*tN,(i-2)*tN+1:(i-1)*tN) = -(2*i-1)/2*I;

A = 1/rstep^2 * A;
% Right hand term
for i = 1:rN-1
    for j = 1:tN
        P(i,j) = i * pressureRT(i*rstep,j*tstep);
    end
end
P(1,:) = P(1,:) + 1/8  * pressureXY(0,0) * ones(1,tN);
%% Solving Linear System by CG
for i = 1:rN-1  
    for j = 1:tN
        % 1 step Jacobian Iteration result as the initial value
        U(i,j) = P(i,j)/A((i-1)*tN+j,(i-1)*tN+j);
    end
end
r0 = P - SparseAx(A,U,rN,tN); p = r0;  % SparseAx is for calculate Ax faster
% begin CG
alpha = 0; beta = 0; k = 0; 
while max(max(abs(r0))) > 0.01  % standard of stopping the iteration
    V = SparseAx(A,p,rN,tN);
    alpha = sum(sum(r0.*r0))/sum(sum(p.*V));
    U = U + alpha * p;
    r1 = r0 - alpha * V;
    beta = sum(sum(r1.*r1))/sum(sum(r0.*r0));
    r0 = r1;
    p = r0 + beta * p;
    k = k + 1;
    % print the residue's infinity norm
    disp(['after ',num2str(k),' steps, the infty norm residue is ',num2str(max(max(abs(r0))))])
end

%% Plot the Deflection 
%(uncommenting the lines including "contourf" to get the contour figure.)
UData = zeros(rN+1,tN+1);
UData(2:rN,2:tN+1) = U; UData(2:rN,1) = U(:,tN); 
UData(1,:) = (sum(U(1,:))/tN + 1/4 * rstep^2 * pressureXY(0,0)) * ones(1,tN+1);

theta = tstep * (0:tN);
r = rstep * (0:rN);
x=r'.*cos(theta);
y=r'.*sin(theta);
surf(x,y,UData);          % get the surface of Deflection
% contourf(x,y,UData,10); % get the contour of Deflection
% axis equal;
%% Plot the Load
%(uncommenting the lines including "contourf" to get the contour figure.)
LdData = zeros(rN+1,tN+1);
for i = 2:rN+1
    for j = 1:tN+1
        LdData(i,j) = pressureRT((i-1)*rstep,(j-1)*tstep);
    end
end
for j = 1:tN+1
    LdData(1,j) = pressureXY(0,0);
end

theta = tstep * (0:tN);
r = rstep * (0:rN);
x=r'.*cos(theta);
y=r'.*sin(theta);
surf(x,y,LdData);          % get the surface of Load
% contourf(x,y,LdData,10); % get the contour of Load
% axis equal;
%% Plot the Deflection and Load on y-axis
Uy = zeros(1,2*rN+1);
Ldy = zeros(1,2*rN+1);
for j = 1:rN
    Uy(rN+1+j) = UData(j+1,tN/4+1);
    Uy(rN+1-j) = UData(j+1,3*tN/4+1);
    Ldy(rN+1+j) = LdData(j+1,tN/4+1);
    Ldy(rN+1-j) = LdData(j+1,3*tN/4+1);
end
Uy(rN+1) = UData(1,1);
Ldy(rN+1) = LdData(1,1);
yaxis = rstep * (-rN:rN);
plot(yaxis,Ldy);   % plot the load on y-axis
hold on;
plot(yaxis,50*Uy); % plot the Deflection on y-axis ( scalar '50' can be changed.)