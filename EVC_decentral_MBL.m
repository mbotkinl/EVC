%EVC decentralized PWL
%Micah Botkin-Levy
%Spring 2018
clc;clear all;

load('EVCscenarioN6.mat') %generated with "EVC_scenario_MBL.m"

%desired states
Sn0=SOCmin;
Kn=FullChargeTime;

%initialize
lambdaGuess=1;
lambda0=ones(K+1,1)*lambdaGuess;
lambda=lambda0;
alpha=.5;
numIteration=100;

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
            
            if cvx_status == "Failed"
                fprintf("Optimization Failed")
                break
            else
                Xn{p,1}(2:K+2,evInd)=xn; %solved state goes in next time slot
                U{p,1}(:,evInd)=un;  %current goes in current time slot
            end
        end
        
        %solve coordinator problem
        cvx_begin quiet
            variable z(S*(K+1),1)
            variable xt(K+1,1)
            minimize (-sum(lambda)*sum(z))
            subject to
                (eye(K+1)-AhatT)*xt==AhatT0*T0+VhatT*w+EhatT*z;
                xt<=Tmax;
                xt>=0;
                z>=0;
                z<=deltaI;
        cvx_end
        
        if cvx_status == "Failed"
            fprintf("Coordinator Optimization Failed")
            break
        else
            Xt(:,p)=xt;
        end


        %grad of lagragian
        gradL=sum(U{p,1},2)+Ghat*w-Fhat*z;
        %update lambda
        lambda_new=lambda+alpha*gradL;
        Lam(:,p)=lambda_new;
        lambda=lambda_new;
    end
%end



