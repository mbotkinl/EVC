clc;clear all;
EV_Charging_Scenario
x0=[s0;T0];
%fix this repeat value to get K+1 ???
Dload=[Dload';Dload(end)];
Tamb=[Tamb;Tamb(end)];
Vload=34.5; %kV
w=zeros(2*(K+1),1);
ind=1;
for i=1:K+1
w(ind)=Dload(i);
w(ind+1)=Tamb(i);
ind=ind+2;
end
%???
M=S;
%create matrices
Ahat=[[zeros(N+1,(K)*(N+1));blkdiagMat(Ad,K)] zeros((K+1)*(N+1),N+1)];
A0hat=[Ad;zeros((N+1)*(K+1)-(N+1),N+1)];
Bhat=blkdiagMat(Bd,K+1);
Vhat=blkdiagMat(Vd,K+1);
Ehat=blkdiagMat(Ed,K+1);
Hhat=blkdiagMat(Ha,K+1);
Ghat=blkdiagMat(Ga,K+1);
Fhat=blkdiagMat(Fa,K+1);
%initialize MPC
%finalT=ceil(1.2*max(FullChargeTime));
%finalT=ceil(5*max(FullChargeTime));
%finalT=Ts*50;
steps= 30;  %round(finalT/Ts);
X=zeros((N+1),steps);
X(:,1)=x0;
xi=x0;
U=zeros(N,steps);
%i=1;
for i=1:steps
fprintf("MPC step %g of %g....\n",i,steps)
%solve optimization from step i to K+i
cvx_solver Gurobi
%cvx_solver SDPT3
%cvx_precision low
cvx_begin quiet
variable u(N*(K+1),1) %control
variable x((N+1)*(K+1),1) %states (time 1 all states first N+1 rows)
%variable w(2*(K+1),1)  %distrubances
variable z(M*(K+1),1) %pw currents
minimize (u'*Rt*u+x'*Qt*x-2*ones(1,(N+1)*(K+1))*Qt*x)
subject to
(eye((K+1)*(N+1))-Ahat)*x==A0hat*xi+Bhat*u+Vhat*w+Ehat*z; %64
0==Hhat*u+Ghat*w+Fhat*z; %65
x<=repmat([ones(N,1);Tmax],K+1,1);
x>=0;
u<=repmat(imax,K+1,1);
u>=repmat(imin,K+1,1);
z>=0;
z<=deltaI;
cvx_end
%     d=x(1:N+1:length(x));
%     MPCfig(d,i)
%     pause
%x0(N+1)
%x(N+1)
%A0hat(N+1,:)*xi+Vhat(N+1,:)*w+Ehat(N+1,:)*z
%     for ii= 1%:N
%         figure
%         subplot(2,1,1)
%         plot(x(ii:N+1:length(x)))
%         subplot(2,1,2)
%         plot(u(ii:N:length(u)))
%     end
%     figure
%     plot(x(N+1:N+1:length(x)))
%     pause
%     close all
%apply u1 and get new xi for i+1
xi=x(1:N+1);
X(:,i+1)=xi;
U(:,i)=u(1:N);
end