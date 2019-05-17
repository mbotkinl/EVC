#generate simulation scenario for Centralized EV Charging Problem
#Micah Botkin-Levy
#4/13/18
#from Mads Almassakhi code

function setupScenario(N;Tmax=100,num_homes=0,Dload_error=0,saveS=false,path=pwd())

    #model parameters
    a   = rand(N,1)*.1 .+ 0.8        # efficiency of Li-ion batts is ~80-90%
    b_options = [40 60 75 100]       # vehicle battery choices (kWh)
    b_prob = [.2 .4 .85]             # distribution of batterys
    r_ind = rand(N,1)
    b_kWh = b_options[1]*ones(Int,N,1)       # vector of battery sizes (kWh)
    b_kWh[r_ind.>b_prob[1]] .= b_options[2]
    b_kWh[r_ind.>b_prob[2]] .= b_options[3]
    b_kWh[r_ind.>b_prob[3]] .= b_options[4]
    b=b_kWh*3.6e6 # battery capacity (MJ)

    #transformer parameters
    xfrmR  = 10e3/3                        # single phase transformer rating kVA
    Vac = 240                              # PEV battery rms voltage --- V [used in PEV kW -> kA conversion]
    Vtf = 8320                             # distr-level transformer rms voltage --- V [used in inelastic kW -> kA conv]
    Ntf   = Vtf/Vac                        # pole-top transformer turns ratio

    # New Transformer Dynamics Paramters scaled from Hines paper
    power_weight = 0.000939/(1e3/3/2.5)^2*1.5/60
    curr_weight = power_weight*Vac^2
    beta = 0.0149*2/60
    alpha = 0.178*5/60

    # transformer parameters
    Ts=3*60 #seconds
    τP = exp(- Ts*beta)
    ρP = 1 - τP
    γP = 1/beta*ρP*curr_weight

    # vehicle energy dynamic parameters
    ηP = round.(Ts*Vac*a./b*1000,digits=4)  # 1/kA, normalized battery sizes (0-1)

    # PWL Parameters
    S=6                              # number of PWL segments
    ItotalMax =  (xfrmR/Vtf)*Ntf*1.8 # kA can overload by 1.8 p.u
    deltaI = ItotalMax/S             # width of PWL segments

    # MPC Paramters
    T1=8
    T2=14-T1
    K1 = round(Int,T1*3600/Ts);        # Initial Prediction and Fixed Horizon (assume K1 instants = 12 hrs)
    K2 = round(Int,T2*3600/Ts);        # Additional time instants past control horizon
    K  = K1+K2;                        # Total horizon (8 PM to 10 AM)

    # Constraint parameters:
    imin = zeros(N,1)                                # enable v2g with imin<0
    b_nom = (b_kWh .- minimum(b_options))/(maximum(b_options)-minimum(b_options))
    shift=2
    i_nom = min.(max.((b_nom .+ rand(Beta(4,3),N,1)/shift .- 1/(2*shift)),0),1)
    imax =round.(((80-28)*i_nom.+28)/1000,digits=6)  # kA, charging with 28-80 A

    # Initial conditions:
    s0 = round.(((.7)*rand(Beta(4, 3),N,1)),digits=4)   # initial states of charge
    t0=Tmax-30                                          # initial temp (~30 C below Tmax)

    #desired states
    Snmin = round.(((1-.75)*rand(Beta(3, 2),N,1).+.75),digits=4)   # Required min final states of charge (~0.75-1)
    FullChargeTime_relative = round.(((1-.75)*rand(Beta(2, 3),N,1).+.75),digits=4)
    Kn = convert(Array{Int,2},round.(K*FullChargeTime_relative))   # time step for final state of charge

    # Disturbance scenario
    peak_demand_house = 4200 #W
    min_demand_house = 3200 #W
    Dload_amplitude=num_homes*peak_demand_house
    Dload_minimum = num_homes*min_demand_house
    dist = [range(0,stop=10,length=Int(round(K/2)));range(-10,stop=-1,length=Int(K-round(K/2)))] # let demand per household be peaking at 8PM and 8 PM
    d = Normal(0,3)
    inelasticDemand = pdf.(d,dist)
    FullinelasticDemand = (inelasticDemand.-minimum(inelasticDemand))/(maximum(inelasticDemand)-minimum(inelasticDemand))
    FullDload=reshape((Dload_amplitude-Dload_minimum)*FullinelasticDemand.+Dload_minimum ,(K,1));    # total non-EV demand (in W)
    iD_pred = round.(FullDload/Vac/1e3,digits=6)    #background predicted demand current (kA)
    noisePerc= Dload_error/Dload_amplitude
    iD_actual = round.(iD_pred+2*noisePerc*iD_pred.*rand(length(iD_pred),1).-iD_pred*noisePerc,digits=6) #background actual demand current (kA)

    Tamb_raw  = round.(Tamb_amplitude*ones(K+1,1).-pdf.(d,range(-8,stop=8,length=K+1))*10,digits=6); # ambient temperature (C)
    Tamb = Tamb_raw.+alpha/beta    # adding constant dynamic from Hines paper

    # penalty matrix
    qrRatio=round.((500-0.005)*rand(Beta(0.5, 0.5),N,1) .+ 0.005,digits=2)
    Ri=1e1*ones(N,K)
    Qi=qrRatio.*Ri/100

    # reduce Q and R for time steps after vehicle can leave
    for n=1:N
        Ri[n,Kn[n]:K].=Ri[n,1]*100
        Qi[n,Kn[n]:K].=Qi[n,1]/10
    end
    QT=zeros(1,K)
    Qsi=[Qi;QT];

    #for slack (not used)
    β=1e3*rand(N,1)

    @assert all(ηP.*Kn.*imax+s0 .>= Snmin) "Some PEVs may not be able to meet SOC min level by desired time!"
    evS=scenarioStruct(N,Ts,K1,K2,K,S,ItotalMax,deltaI,Tmax,imin,imax,a,b_kWh,
                       ηP,τP,ρP,γP,s0,t0,Snmin,Kn,iD_pred,iD_actual,Tamb,Tamb_raw,Qsi,Ri,β)

    if saveS==true
        save(path*"EVCscenarioN$(N).jld2","evScenario",evS)
    end

    return evS
end

function setupHubScenario(H,Nh;Tmax=.393,Dload_amplitude=0,saveS=false,path=pwd())
    #model parameters
    a   = rand(1,H)*.1 .+ 0.8                  # efficiency of collect hubs
    b_options = [100 200 600]                  # battery size options
    b_prob = [0 .4 .8; 1 0 .7; 1 0 1;1 0 .5]   # distributions in each hub - hard coded for now
    r_ind = rand(Nh,1)
    b_kWh = b_options[1]*ones(Int,Nh,H)        # battery capacity (kWh)
    for h=1:H
        b_kWh[(r_ind.>b_prob[h,1])[:],h] .= b_options[1]
        b_kWh[(r_ind.>b_prob[h,2])[:],h] .= b_options[2]
        b_kWh[(r_ind.>b_prob[h,3])[:],h] .= b_options[3]
    end
    b=b_kWh*3.6e6 # battery capacity (MJ)
    imax = round.(((1-.2)*rand(Beta(3, 2),Nh,H).+.2),digits=4) # kA, charging with 200-1000 A

    # transformer scneario
    xfrmR  = 120e3/3                        # single phase transformer rating kVA
    Vac = 480                               # PEV battery rms voltage --- V [used in PEV kW -> kA conversion]
    Vtf = 13.2e3                            # distr-level transformer rms voltage --- V [used in inelastic kW -> kA conv]
    Ntf   = Vtf/Vac                         # pole-top transformer turns ratio

    # Discretization parameters scaled from Hines paper
    Ts=180.0
    power_weight = 0.000939/(1e3/3*4)^2*1.5/60
    curr_weight = power_weight*Vac^2
    beta = 0.0149*2/60
    alpha = 0.178*5/60

    # transformer dynamic parameters
    Ts=3*60 #seconds
    τP = exp(- Ts*beta)
    ρP = 1 - τP
    γP = 1/beta*ρP*curr_weight

    ηP=a*(Ts/3600)*Vac/1e3   # 1/kA, normalized battery sizes (0-1)

    ## MPC Paramters
    T1=12 #hours
    T2=14-T1
    K1 = round(Int,T1*3600/Ts);            # Initial Prediction and Fixed Horizon (assume K1 instants = 12 hrs)
    K2 = round(Int,T2*3600/Ts);             # Additional time instants past control horizon
    K  = K1+K2;                        # Total horizon (8 PM to 10 AM)

    # PWL Parameters:
    S=15
    ItotalMax =  (xfrmR/Vtf)*Ntf*1.8 # kA can overload by 1.8 p.u
    deltaI = ItotalMax/S             # width of PWL segments

    #system information
    t0=Tmax-30
    Tamb_amplitude = 18 #C

    # Disturbance scenario
    Dload_error=0
    num_homes=1
    peak_demand_house = 30e6  #W
    min_demand_house =25e6 #W
    Dload_amplitude=num_homes*peak_demand_house #W
    Dload_minimum = num_homes*min_demand_house
    dist = [range(-1,stop=10,length=Int(round(K/2)));range(-10,stop=-1,length=Int(K-round(K/2)))]
    d = Normal(0,3)
    inelasticDemand = pdf.(d,dist)
    FullinelasticDemand = (inelasticDemand.-minimum(inelasticDemand))/(maximum(inelasticDemand)-minimum(inelasticDemand))
    FullDload=reshape((Dload_amplitude-Dload_minimum)*FullinelasticDemand.+Dload_minimum ,(K,1));    # total non-EV demand (in W)
    iD_pred = round.(FullDload/Vac/1e3,digits=6)    #background demand current (kA)
    if Dload_error>0
        noisePerc= Dload_error/Dload_amplitude
        iD_actual = round.(iD_pred+2*noisePerc*iD_pred.*rand(length(iD_pred),1).-iD_pred*noisePerc,digits=6)
    else
        iD_actual=iD_pred
    end
    Tamb_raw  = round.(Tamb_amplitude*ones(K+1,1).-pdf.(d,range(-8,stop=8,length=K+1))*10,digits=6);
    Tamb = Tamb_raw.+alpha/beta

    # penalty matrix
    qrRatio=round.((200-0.1)*rand(Beta(0.6, 0.6),1,H) .+ 0.1,digits=2)
    orRatio=10
    Rh=1e1*ones(1,H)
    Qh=qrRatio.*Rh/10
    Oh=orRatio.*Rh
    #action happens immediately afte interval before ends

    #hub information
    arriveInd = rand(1:7,H) # last hour of arrival
    departInd = rand(8:11,H) # first hour of departure
    K_arrive_pred=zeros(Nh,H) # predicted time step of arrival
    K_depart_pred=zeros(Nh,H) # predicted time step of departure
    for h=1:H
        arriveLast=round(Int,arriveInd[h]*3600/Ts) # time step of last arrival
        departFirst=round(Int,departInd[h]*3600/Ts) # time step of first departure
        # these beta functions change the distribution of vehicles arriving and leaving
        K_arrive_pred[:,h]=rand(BetaBinomial(Int(arriveLast*1.5),rand(1:8),rand(1:8)),Nh,1) .- Int(arriveLast/2)
        K_depart_pred[:,h]=rand(BetaBinomial(K-departFirst,rand(1:8),rand(1:8)),Nh,1) .+ departFirst
    end

    # for now there is no forecast error with arrival and departure times
    K_arrive_actual=K_arrive_pred
    K_depart_actual=K_depart_pred
    Sn_depart_min=round.(1 .- 0.20*rand(Nh,H),digits=4) # minimum departure energy: need 80-100%
    Sn_arrive_pred=round.((0.4-0.1)*rand(Nh,H) .+ 0.1 ,digits=4) #arrival energy  10-40%
    Sn_arrive_actual=Sn_arrive_pred
    EVcap=round.(b./3.6e6/1e3,digits=6) #MWh
    e0=zeros(H) # initial energy in each hub (from vehicle already arrived)
    for h=1:H
        ind= findall(x->x==0,K_arrive_actual[:,h])
        if isempty(ind)
            e0[h]=0
        else
            e0[h]=sum(Sn_arrive_actual[n,h]*EVcap[n,h] for n in ind)
        end
    end

    #prepare predicted values for optimization
    eMax=zeros(K,H)             # maximum energy capacity
    uMax=zeros(K,H)             # maximum chargu=ing current
    eDepart_min=zeros(K,H)      # minimum energy departing
    eArrive_pred=zeros(K,H)     # predicted energy arriving
    eArrive_actual=zeros(K,H)   # actual energy arriving
    slackMax=zeros(K,H)         # maximum "extra" energy to fill between minimum energy departing and energy capacity departing
    for k =1:K
        for h=1:H # find a way to do this without another loop???
            depart=[n for n=1:Nh if k==K_depart_pred[n,h]]
            if length(depart)!=0
                eDepart_min[k,h] = sum(Sn_depart_min[n,h]*EVcap[n,h] for n in depart)
                slackMax[k,h] = sum(EVcap[n,h]-Sn_depart_min[n,h]*EVcap[n,h] for n in depart)
            end
            arrive=[n for n=1:Nh if k==K_arrive_pred[n,h]]
            if length(arrive)!=0
                eArrive_pred[k,h] = sum(Sn_arrive_pred[n,h]*EVcap[n,h] for n in arrive)
                eArrive_actual[k,h] = sum(Sn_arrive_actual[n,h]*EVcap[n,h] for n in arrive)
            end
            parked=[n for n=1:Nh  if K_arrive_pred[n,h]<=k<K_depart_pred[n,h]]
            if length(parked)!=0
                #ηH[k]=mean(eta[n] for n in parked)
                uMax[k,h]=sum(imax[n,h] for n in parked)
                eMax[k,h]=sum(EVcap[n,h] for n in parked)
            end
        end
    end

    hubS=scenarioHubStruct(Nh,H,Ts,K1,K2,K,S,ItotalMax,deltaI,Tmax,ηP,τP,ρP,γP,e0,t0,iD_pred,iD_actual,Tamb,Tamb_raw,Qh,Rh,Oh,
                            Sn_depart_min,Sn_arrive_actual,Sn_arrive_pred,K_arrive_pred,K_depart_pred,K_arrive_actual,
                            K_depart_actual,EVcap,eMax,uMax,eDepart_min,eArrive_pred,eArrive_actual,slackMax)

    if saveS==true
            save(path*"HubscenarioH$(H).jld2","hubS",hubS)
    end

    return hubS
end
