%generate simulation scenario for Centralized EV Charging Problem
%Micah Botkin-Levy
%3/26/18
%from Mads Almassakhi code
clc;clear all;
N   = 6;
%% model parameters 1:
a = rand(N,1)*.1 + 0.8; % efficiency of Li-ion batts is ~80-90%
b = 6*rand(N,1)+12;     % battery capacity (12-18 kWh)
eta = a./b;             % normalized battery sizes (0-1)
Rw = 1.0;               % Ohm
m = 13;                 % transformer mass in kg           -------> CAUTION: This is way too low!!! Physics not right!!!
%m=555;
C = 0.45*m;             % heat cap. thermal mass ----- spec. heat cap. of {carbon steel, iron, veg. oil} = {0.49, 0.45, 1.67} kJ/(kg*K)
R = 1070e-4/(5*(m/7870)^(2/3)); % heat outflow resistance K/kW : R = 0.1073 (K*m^2/W)/(A_s), rule of thumb calculation
Ts = 155;                 % sampling time in seconds
%Ts=.155;
tau = 1 - Ts/(R*C);     % temp time constant: 1 - 1/RC
rho = 1 - tau;          % ambient-to-temp param: 1/RC
gamma = Ts*Rw/C;           % load-to-temp param:
Vtf = 8320;
%% Model parameters 2:
% a   = rand(N,1)*.1 + 0.8;               % efficiency of Li-ion batts is ~80-90%
% b   = (6*rand(N,1)+12)*3.6e6;           % battery capacity (12-18 kWh = 43.3-64.8 MJ)
% m   = 2000;                             % transformer mass in kg
% C   = 450*m;                            % heat cap. thermal mass J/K ----- spec. heat cap. of C = {carbon steel, iron, veg. oil} = {490, 450, 1670} J/(kg*K)
% Rh   = 1070e-4/(35*5*(m/7870)^(2/3));    % heat outflow resistance K/W : R = 0.1073 (K*m^2/W)/(A_s), rule of thumb calculation
% Rw  = 1;                              % coil winding resistance --- ohms:
% Vac = 240;                              % PEV battery rms voltage --- V [used in PEV kW -> kA conversion]
% Vtf = 8320;                             % distr-level transformer rms voltage --- V [used in inelastic kW -> kA conv]
% Ntf   = Vtf/Vac;                          % pole-top transformer turns ratio
% %% Discretization parameters:
% Ts = Rh*C/9;              % s, sampling time in seconds
% eta = Ts*Vac*a./b;        % 1/A, normalized battery sizes (0-1)
% tau = 1 - Ts/(Rh*C);       % no units, temp time constant: 1 - 1/RC
% rho = 1 - tau;            % no units, ambient-to-temp param: 1/RC
% gamma = Ts*Rw/(C*Ntf);      % K/W, ohmic losses-to-temp parameter
%% PWL Parameters:
S = 3;
ItotalMax = 20;        % CAUTION  ---> Imax gives upper limit on total current input on Transfomer and if picked too low will cause infeasible.
deltaI = ItotalMax/S;
%% MPC Paramters
%K1 = round(12*3600/Ts);            % Initial Prediction and Fixed Horizon (assume K1 instants = 12 hrs)
%K2 = round(2*3600/Ts);             % Additional time instants past control horizon
%K  = K1+K2;                        % Total horizon (8 PM to 10 AM)Qs  = 10;              % Stage and terminal penalty on charge difference with respect to 1 (states s)
K=99;
%% Constraint parameters:
Tmax = 393;                             % Short-term over-loading --> 120 C = 393 Kelvin
imin = zeros(N,1);                      % A, q_min < 0 if V2G is allowed
imax = (10 + 16*rand(N,1));             % A, charging with 10-24 A
SOCmin = 1 - 0.20*rand(N,1);            % Required min final states of charge (~0.80-1)
FullChargeTime_relative = .25*rand(N,1)+.75;
FullChargeTime = round(K*FullChargeTime_relative);
% Initial conditions:
s0 = 0.25*rand(N,1);      % initial states of charge (0 - 0.20)
T0 = 368;                 % initial temp (~65 K below Tmax)
% Disturbances
Dload_amplitude = 2;  % base-demand factor
Tamb_amplitude  = 303;   % assume hot night in summer (30 C)
%% constraint matrices
Et=gamma*deltaI*[1:2:(2*S-1)];
% Augmented model (see notes, for simulation only) these are from 43-46
% x[k+1] = Ad*x[k] + Bd*u[k] + Vd*w[k] + Ed*z[k]
%      0 =           Ha*u[k] + Ga*w[k] + Ha*z[k]
Ad = diag([ones(N,1);tau]);
Bd = [diag(eta);zeros(1,N)];
Vd = [zeros(N,2); 0, rho];
Ed = [zeros(N,S);Et];
Fa = ones(1,S);
Ga = [-1 0];
Ha = -ones(1,N);
%create matrices for compact central problem
Ahat=[[zeros(N+1,(K)*(N+1));blkdiagMat(Ad,K)] zeros((K+1)*(N+1),N+1)];
A0hat=[Ad;zeros((N+1)*(K+1)-(N+1),N+1)];
Bhat=blkdiagMat(Bd,K+1);
Vhat=blkdiagMat(Vd,K+1);
Ehat=blkdiagMat(Ed,K+1);
Hhat=blkdiagMat(Ha,K+1);
Ghat=blkdiagMat(Ga,K+1);
Fhat=blkdiagMat(Fa,K+1);
%create matricies for dual/decentral
Ahats=diag(ones(K,1),-1);
Ahats0=[1;zeros(K,1)];
Bhats=cell(N,1);
for i=1:N
Bhats{i,1}=diag(ones(K+1,1)*eta(i));
end
AhatT=diag(ones(K,1)*tau,-1);
AhatT0=[tau;zeros(K,1)];
VhatT=blkdiagMat([0 rho],K+1);
EhatT=blkdiagMat(Et,K+1);
%% Disturbance scenario:
%FullinelasticDemand = [normpdf(0,linspace(0,8,round((K1-1)/2)),3) normpdf(0,linspace(-8,0,round(K1/2)),3)]; % let demand per household be peaking at 8PM and 8 PM with nadir inbetween
%FullinelasticDemand = 100*(200*(FullinelasticDemand-min(FullinelasticDemand))/range(FullinelasticDemand) + 600)/1000; % total non-EV demand (in kW) = N/PEVpenetration*inelasticDemandperHouse
%FullinelasticDemand = [FullinelasticDemand'; FullinelasticDemand(end)*ones(K2+1,1)];
inelasticDemand = [normpdf(0,linspace(0,8,round(K/2)),3) normpdf(0,linspace(-8,0,round(K/2)),3)]; % let demand per household be peaking at 8PM and 8 PM with nadir inbetween
FullinelasticDemand = 100*(200*(inelasticDemand-min(inelasticDemand))/range(inelasticDemand) + 600)/1000; % total non-EV demand (in kW) = N/PEVpenetration*inelasticDemandperHouse
FullDload   = Dload_amplitude*FullinelasticDemand;    % peaks during mid-day
iD=FullDload/Vtf; %background demand current
FullTamb    = Tamb_amplitude*ones(K+1,1); %normpdf(0,linspace(-10,10,max(K,kmax)),3)';   % exogenous peaks during mid-day          % OVER-NIGHT CHARGING: TIMES -1?
ndisturbs=2;
w = zeros((K+1)*ndisturbs,1);
for i=1:K+1
w((i-1)*ndisturbs+1:i*ndisturbs,1)  = [iD(i); Tamb_amplitude];
end
%% penalty matrix
Qs  = 10;              % Stage and terminal penalty on charge difference with respect to 1 (states s)
QT  = 0;               % PENALTY ON TEMPERATURE DEVIATION (W.R.T 0)
R   = .1;              % Stage and terminal penalty on local power flow (inputs q)
QsKi  = 1;             % Stage and terminal penalty on charge difference with respect to 1 (states s), for k >= Ki
RKi   = 10;            % Stage and terminal penalty on local power flow (inputs q), for k >= Ki
Qsi = ones(N,1);
Ri = ones(N,1);
QsiKi = ones(N,1);
RiKi = ones(N,1);
for i=1:N
Qsi(i)   = Qs   * (10*rand(1)+.01);
Ri(i)    = R    * (5*rand(1)+.1);
QsiKi(i) = QsKi * (rand(1)+.1);
RiKi(i)  = RKi  * (rand(1)+.1);
end
Rt = []; Qt = [];
for k=1:K+1
%for k=1:K
for i=1:N
if k<FullChargeTime(i)
Rt = blkdiag(Rt,  Ri(i));
Qt = blkdiag(Qt, Qsi(i));
else
Rt = blkdiag(Rt,  RiKi(i));
Qt = blkdiag(Qt, QsiKi(i));
end
end
Qt = blkdiag(Qt, QT);    % Temperature penalty
end
Rt = sparse(Rt); Qt = sparse(Qt);
%% save
if( any(eta.*K.*FullChargeTime_relative.*imax+s0 < SOCmin) )
%if( any(eta.*K1.*FullChargeTime_relative.*imax+s0 < SOCmin) )
disp('Some PEVs may not be able to meet SOC min level by desired time!');
end
name=sprintf("EVCscenarioN%d.mat",N);
save(name)