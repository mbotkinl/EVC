%EVC decentralized PWL
%Micah Botkin-Levy
%Spring 2018
clc;clearvars;

N=20;
Testnum=1;
testFolder=sprintf('N%d_T%d',N,Testnum);
scenarioFile=sprintf('EVCscenarioN%d.mat',N);
load(fullfile(testFolder,scenarioFile)) %generated with EVC_scenario_MBL

%desired states
Sn0=SOCmin;
Kn=FullChargeTime;

%initialize
lambdaGuess=10;
lambda0=ones(K+1,1)*lambdaGuess;
lambda=lambda0;
alpha=.5;
numIteration=500;
steps=K+1;

%for step=1:steps
    step=1;
    Lam=zeros((K+1),numIteration); %(rows are time, columns are iteration)
    Lam(:,1)=lambda0;
    Xt=zeros((K+1),numIteration); %(rows are time, columns are iteration)
    Xn=cell(numIteration,1); %each cell is a new iteration
    for i=1:numIteration
        Xn{i,1}=zeros(K+2,N); %each cell has column for each EV, rows are time
        Xn{i,1}(1,:)=s0';
    end
    
    U=cell(numIteration,1); %each cell is a new iteration
    for i=1:numIteration
        U{i,1}=zeros(K+1,N); %each cell has column for each EV, rows are time
    end
    
    %p=1;
    for p=1:numIteration
        

        fprintf("iteration step %g of %g....\n",p,numIteration)
        %solve N subproblems
        
        
        %xtemp=zeros(K+2,N);
        %utemp=zeros(K+1,N);
        %parfor evInd=1:N
        
        for evInd=1:N
            %move this elsewhere after (but  need to change as i-->K)
            Qhatn=eye(K+1)*Qsi(evInd);
            Rhatn=eye(K+1)*Ri(evInd);
            
            cvx_solver Gurobi
            
            cvx_begin quiet
                target=zeros((K+1),1);
                target(max(1,Kn(evInd)-(step-1)*Ts):length(target),1)=Sn0(evInd); %fix Ts for time loop???
                variable xn(K+1,1)
                variable un(K+1,1)
                minimize (un'*Rhatn*un+xn'*Qhatn*xn-2*ones(1,(K+1))*Qhatn*xn+lambda'*un)
                subject to
                    (eye(K+1)-Ahats)*xn==Ahats0*s0(evInd)+Bhats{evInd,1}*un;
                    xn<=1;
                    xn>=target;
                    un<=imax(evInd);
                    un>=imin(evInd);
            cvx_end
            
            if cvx_status ~= "Solved"
                fprintf("Optimization Failed")
                %break
            else
                Xn{p,1}(2:K+2,evInd)=xn; %solved state goes in next time slot
                %xtemp(:,evInd)=xn;
                U{p,1}(:,evInd)=un;  %current goes in current time slot
                %utemp(:,evInd)=un;
            end
        end
        
        
        %solve coordinator problem
        cvx_begin quiet
            variable z(S*(K+1),1)
            variable xt(K+1,1)
            %minimize (-sum(lambda)*sum(z))
            %minimize (-sum(lambda(1:K+1,1))*sum(z(1:S,1)))
            minimize (-lambda'*Fhat*z)
            subject to
                (eye(K+1)-AhatT)*xt==AhatT0*T0+VhatT*w+EhatT*z;
                xt<=Tmax;
                xt>=0;
                z>=0;
                z<=deltaI;
        cvx_end
        
        if cvx_status ~= "Solved"
            fprintf("Coordinator Optimization Failed")
            %break
        else
            Xt(:,p)=xt;
        end


        %grad of lagragian
        gradL=sum(U{p,1},2)+Ghat*w-Fhat*z;
        
        %update lambda
        alpha_p = alpha/ceil(p/10);
        
        %alpha_p=alpha;
        lambda_new=lambda+alpha_p*gradL;
        Lam(:,p)=lambda_new;
        lambda=lambda_new;
        
        
    end
%end


figure; hold on;
for ii =1
    plot(Lam(ii,:))
end
% plotName='Lam1';
%print(fullfile(testFolder,plotName),'-dpng','-r0')

figure;
lgd=strings(1,N);
for ii =1:N
    subplot(3,1,1)
    plot(Xn{numIteration,1}(:,ii));
    hold on;
    ylabel("SOC")
    xlim([1 steps])
    subplot(3,1,2)
    plot(U{numIteration,1}(:,ii));
    hold on;
    ylabel("Current")
    xlim([1 steps])
    lgd(ii)=sprintf("EV%d",ii);
end
legend(lgd)

subplot(3,1,3)
plot(Xt(:,numIteration))
ylabel("XFRM Temp (K)")
xlim([1 steps])

plotName='Decentral1';
%print(fullfile(testFolder,plotName),'-dpng','-r0')

