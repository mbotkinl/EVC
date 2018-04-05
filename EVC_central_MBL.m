%EVC centralized compact PWL
%Micah Botkin-Levy
%Spring 2018
clc;clearvars;

N=10;
Testnum=1;
testFolder=sprintf('N%d_T%d',N,Testnum);
scenarioFile=sprintf('EVCscenarioN%d.mat',N);
load(fullfile(testFolder,scenarioFile))  %generated with EVC_scenario_MBL

%initialize
steps= K+1;  %round(finalT/Ts);
X=zeros((N+1),steps);
X(:,1)=x0;
xi=x0;
U=zeros(N,steps);
cvx_solver Gurobi

i=1;


%for i=1:steps
    fprintf("MPC step %g of %g....\n",i,steps)
    %solve optimization from step i to K+i
    %cvx_solver SDPT3
    %cvx_precision low
    %cvx_solver_settings('NumericFocus',2)
    %cvx_solver_settings('ResultFile','test.lp')

    %gurobi_cl C:/Users/micah/Documents/uvm/Research/EVC/test.lp
    
    cvx_begin %quiet
        target=zeros((N+1)*(K+1),1); 
        for ii=1:N
           cur=Kn(ii)-(i-1);
           ind=max(0,(cur-1)*(N+1))+ii:N+1:length(target);
           target(ind)=Sn(ii);
        end
        variable u(N*(K+1),1) %control
        variable x((N+1)*(K+1),1) %states (time 1 all states first N+1 rows)
        variable z(S*(K+1),1) %pw currents
        minimize (u'*R*u+x'*Q*x-2*ones(1,(N+1)*(K+1))*Q*x)
        subject to
            (eye((K+1)*(N+1))-Ahat)*x==A0hat*xi+Bhat*u+Vhat*w+Ehat*z; %64
            0==Hhat*u+Ghat*w+Fhat*z; %65
            x<=repmat([ones(N,1);Tmax],K+1,1);
            x>=target;
            u<=repmat(imax,K+1,1);
            u>=repmat(imin,K+1,1);
            z>=0;
            z<=deltaI;
    cvx_end
    %if cvx_status ~= "Solved"
    if cvx_status == "Failed"
        fprintf("Optimization Failed")
        %break
    end
    
    %mpc plot
    %     d=x(1:N+1:length(x));
    %     MPCfig(d,i)
    %     pause
    
    %check iterative step
    %x0(N+1)
    %x(N+1)
    %A0hat(N+1,:)*xi+Vhat(N+1,:)*w+Ehat(N+1,:)*z
    
    %one time step plots
    figure
    lgd=strings(1,N);
    for ii= 1:N
        subplot(3,1,1)
        plot(x(ii:N+1:length(x)))
        ylabel("SOC")
        xlim([1 steps])
        hold on
        subplot(3,1,2)
        plot(u(ii:N:length(u)))
        ylabel("Current")
        xlim([1 steps])
        hold on
        lgd(ii)=sprintf("EV%d",ii);
    end
    subplot(3,1,3)
    plot(x(N+1:N+1:length(x)))
    ylabel("XFRM Temp (K)")
    xlim([1 steps])
    
    
    
    
    
    plotName='Central1';
    print(fullfile(testFolder,plotName),'-dpng','-r0')
    
    
    %apply u1 and get new xi for i+1
    xi=x(1:N+1);s
    X(:,i+1)=xi;
    U(:,i)=u(1:N);
%end




%for MPC format
% figure
% lgd=strings(1,N);
% for ii =1:N
%     subplot(3,1,1)
%     plot(X(ii,:));
%     ylabel("SOC")
%     xlim([1 steps])
%     hold on
%     subplot(3,1,2)
%     plot(U(ii,:));
%     ylabel("Current")
%     xlim([1 steps])
%     hold on
%     lgd(ii)=sprintf("EV%d",ii);
% end
% legend(lgd)
% 
% subplot(3,1,3)
% plot(X(N+1,:))
% ylabel("XFRM Temp (K)")
% xlim([1 steps])
% 
% plotName='Central1';
%print(fullfile(testFolder,plotName),'-dpng','-r0')


