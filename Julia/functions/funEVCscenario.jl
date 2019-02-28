#generate simulation scenario for Centralized EV Charging Problem
#Micah Botkin-Levy
#4/13/18
#from Mads Almassakhi code

function setupScenario(N;Tmax=100,num_homes=0,Dload_error=0,saveS=false,path=pwd())

    #model parameters
    a   = rand(N,1)*.1 .+ 0.8               # efficiency of Li-ion batts is ~80-90%
    b_high=100 #kWh
    b_low = 20 #kWh
    b_kWh=(b_high-b_low)*rand(Beta(1.5, 1.8),N,1).+b_low  # battery capacity (20-100 kWh)
    #histogram(b_kWh,nbins=40)

    # t=.5*rand(Beta(2, 4),N*1000,1)+0.5*rand(Beta(4, 1.2),N*1000,1)
    # t =.5*rand(Normal(38, 10),N*1000,1)+rand(Normal(80, 2),N*1000,1)
    # # use mixture models?
    # histogram(t,nbins=40)


    b=b_kWh*3.6e6 # battery capacity (MJ)
    # a   = 0.8 *ones(N,1)              # efficiency of Li-ion batts is 80%
    # b   = 12*3.6e6                    # battery capacity (12kWh = 43.3-)

    xfrmR  = 10e3/3                          # single phase transformer rating kVA
    # m   = 5e3                             # transformer mass in kg
    # C   = 450*m                            # heat cap. thermal mass J/K ----- spec. heat cap. of C = {carbon steel, iron, veg. oil} = {490, 450, 1670} J/(kg*K)
    # Rh   = 1070e-4/(35*5*(m/7870)^(2/3))   # heat outflow resistance K/W : R = 0.1073 (K*m^2/W)/(A_s), rule of thumb calculation
    # Rw  = 0.01                                # coil winding resistance --- ohms:
    Vac = 240                              # PEV battery rms voltage --- V [used in PEV kW -> kA conversion]
    Vtf = 8320                             # distr-level transformer rms voltage --- V [used in inelastic kW -> kA conv]
    Ntf   = Vtf/Vac                        # pole-top transformer turns ratio

    # Discretization parameters:
    # Ts = Rh*C/9              # s, sampling time in seconds
    # #Ts=180
    # ηP = round.(Ts*Vac*a./b*1000,digits=4)  # 1/kA, normalized battery sizes (0-1)
    # #ηP = Ts*Vac*a./b  # 1/A, normalized battery sizes (0-1)
    #
    # #τP = 1 - Ts/(Rh*C)      # no units, temp time constant: 1 - 1/RC
    # τP = exp(- Ts/(Rh*C))
    #
    # ρP = 1 - τP            # no units, ambient-to-temp param: 1/RC
    # # γP = Ts*Rw/(C*Ntf)*1000^2    # K/kW, ohmic losses-to-temp parameter
    # #γP = Ts*Rw/(C*Ntf)    # K/W, ohmic losses-to-temp parameter
    #
    # #γP = Rh*Rw/(Ntf)*ρP*1000^2    # K/kW, ohmic losses-to-temp parameter
    # γP = Rh*Rw/(Ntf)*ρP*1000^2/1000    # kK/kW, ohmic losses-to-temp parameter
    #     Ts*Rw/(C*Ntf)*1000^2/1000


    # New Discretization Paramters
    power_weight = 0.000939/(1e3/3/2.5)^2*1.5/60
    curr_weight = power_weight*Vac^2
    beta = 0.0149*2/60
    alpha = 0.178*5/60

    Ts=3*60 #seconds
    τP = exp(- Ts*beta)
    ρP = 1 - τP
    γP = 1/beta*ρP*curr_weight
    ηP = round.(Ts*Vac*a./b*1000,digits=4)  # 1/kA, normalized battery sizes (0-1)

    # PWL Parameters:
    #S = 3;
    S=10
    #ItotalMax = 20;        % CAUTION  ---> Imax gives upper limit on total current input on Transfomer and if picked too low will cause infeasible.
    ItotalMax =  (xfrmR/Vtf)*Ntf*1.8 #kA can overload by 1.5 p.u
    #ItotalMax = 4000  #A
    deltaI = ItotalMax/S

    ## MPC Paramters
    T1=10
    T2=14-T1
    K1 = round(Int,T1*3600/Ts);            # Initial Prediction and Fixed Horizon (assume K1 instants = 12 hrs)
    K2 = round(Int,T2*3600/Ts);             # Additional time instants past control horizon
    K  = K1+K2;                        # Total horizon (8 PM to 10 AM)Qs  = 10;

    # Constraint parameters:
    #Tmax = 393                             # Short-term over-loading --> 120 C = 393 Kelvin
    imin = zeros(N,1)                      # A, q_min < 0 if V2G is allowed
    #imax = round.(((80-10)*rand(Beta(3, 6),N,1).+10)/1000,digits=6)  # kA, charging with 10-80 A
    b_nom = (b_kWh.-b_low)/(b_high-b_low)
    i_nom = min.(max.((b_nom .+ rand(Beta(3,6),N,1)/5),0),1)
    # i_nom = min.(max.((b_nom .+ rand(Normal(0,1),N*1000,1)/10),0),1)
    imax =round.(((80-10)*i_nom.+10)/1000,digits=6)
    # histogram(b_nom,nbins=40)
    # histogram(i_nom,nbins=40)
    #histogram(imax,nbins=40)
    #cor(b_kWh,imax)

    #imax=10/1000*ones(N,1)

    # Initial conditions:
    #s0=0.6*ones(N,1)
    s0 = round.(((.7)*rand(Beta(4, 3),N,1)),digits=4)   # initial states of charge
    #histogram(s0,nbins=40)
    #t0 = 320/1000                 # initial temp (~65 K below Tmax) 368K
    t0=40

    #desired states
    Snmin = round.(((1-.75)*rand(Beta(3, 2),N,1).+.75),digits=4)   # Required min final states of charge (~0.80-1)
    FullChargeTime_relative = round.(((1-.75)*rand(Beta(2, 3),N,1).+.75),digits=4)
    #histogram(FullChargeTime_relative,nbins=40)

    #Kn=ones(N,1)
    #Snmin=ones(N,1)
    Kn = convert(Array{Int,2},round.(K*FullChargeTime_relative))

    # Disturbances
    #Dload_amplitude = 2;  # base-demand factor
    #Dload_amplitude = 85 #kWatts?
    #Dload_amplitude = 75000 #Watts?
    #Dload_amplitude = 0
    #Tamb_amplitude  = 313/1000   # assume hot night in summer (30 C) 303K
    Tamb_amplitude = 30 #C

    # Disturbance scenario:
    #num_homes= 1000
    peak_demand_house = 4000 #W
    min_demand_house = 800 #W
    Dload_amplitude=num_homes*peak_demand_house #W
    Dload_minimum = num_homes*min_demand_house

    #dist = [range(0,stop=8,length=Int(round(K/2)));range(8,stop=0,length=Int(K-round(K/2)))] # let demand per household be peaking at 8PM and 8 PM
    dist = [range(-1,stop=10,length=Int(round(K/2)));range(-10,stop=-1,length=Int(K-round(K/2)))] # let demand per household be peaking at 8PM and 8 PM
    d = Normal(0,3)
    inelasticDemand = pdf.(d,dist)
    FullinelasticDemand = (inelasticDemand.-minimum(inelasticDemand))/(maximum(inelasticDemand)-minimum(inelasticDemand))

    #FullinelasticDemand = [FullinelasticDemand; FullinelasticDemand[length(FullinelasticDemand)]*ones(K2+1,1)]
    FullDload=reshape(Dload_amplitude*FullinelasticDemand.+Dload_minimum,(K,1));    # total non-EV demand (in W)
    iD_pred = round.(FullDload/Vac/1e3,digits=6)    #background demand current (kA)
    noisePerc= Dload_error/Dload_amplitude
    iD_actual = round.(iD_pred+2*noisePerc*iD_pred.*rand(length(iD_pred),1).-iD_pred*noisePerc,digits=6)
    Tamb_raw  = round.(Tamb_amplitude*ones(K+1,1).-pdf.(d,range(-10,stop=10,length=K+1))*20,digits=6);             #normpdf(0,linspace(-10,10,max(K,kmax)),3)';   # exogenous peaks during mid-day          % OVER-NIGHT CHARGING: TIMES -1?

    Tamb = Tamb_raw.+alpha/beta

    # stT1=Time(20,0)
    # endT1=Time(23,59)
    # stT2=Time(0,0)
    # endT2=Time(10,0)
    # Xlabels=vcat(collect(stT1:Dates.Second(round(Ts)):endT1),collect(stT2:Dates.Second(round(Ts)):endT2))
    # xticks=(1:40:K,Dates.format.(Xlabels[1:40:K],"HH:MM"))
    # plot(iD_pred*240/num_homes,xticks=xticks)
    # plot(Tamb_raw,xticks=xticks)


    # penalty matrix new (need to fix for k>Ki)
    # Ru   = 0.1*1000^2              # Stage and terminal penalty on local power flow (inputs u)
    # Qs  = 10;               # Stage and terminal penalty on charge difference with respect to 1 (states s)


    Qs=1;
    Ru=Qs*100;
    QT  = 0;                # PENALTY ON TEMPERATURE DEVIATION (W.R.T 0)
    Ri=round.(Ru*(5*rand(N,1).+.1),digits=6);
    Qi=round.(Qs*(10*rand(N,1).+.01),digits=6);
    Qsi=[Qi;QT];

    #for slack
    β=1e3*rand(N,1)

    #move this into struct???
    @assert all(ηP.*Kn.*imax+s0 .>= Snmin) "Some PEVs may not be able to meet SOC min level by desired time!"
    # i=1
    # ηP[i].*Kn[i].*imax[i]+s0[i] .>= Snmin[i]

    evS=scenarioStruct(N,Ts,K1,K2,K,S,ItotalMax,deltaI,Tmax,imin,imax,
                       ηP,τP,ρP,γP,s0,t0,Snmin,Kn,iD_pred,iD_actual,Tamb,Tamb_raw,Qsi,Ri,β)

    if saveS==true
        save(path*"EVCscenarioN$(N).jld2","evScenario",evS)
    end

    return evS
end

function setupHubScenario(H,Nh;Tmax=.393,Dload_amplitude=0,saveS=false,path=pwd())
    #model parameters
    a   = rand(1,H)*.1 .+ 0.8               # efficiency of Li-ion batts is ~80-90%
    b   = (6*rand(Nh,H).+12)*3.6e6           # battery capacity (12-18 kWh = 43.3-64.8 MJ)
    imax = (10 .+ 16*rand(Nh,H))/1000             # kA, charging with 10-24 A

    m   = 2000                             # transformer mass in kg
    C   = 450*m                            # heat cap. thermal mass J/K ----- spec. heat cap. of C = {carbon steel, iron, veg. oil} = {490, 450, 1670} J/(kg*K)
    Rh   = 1070e-4/(35*5*(m/7870)^(2/3))   # heat outflow resistance K/W : R = 0.1073 (K*m^2/W)/(A_s), rule of thumb calculation
    Rw  = 1                                # coil winding resistance --- ohms:
    Vac = 240                              # PEV battery rms voltage --- V [used in PEV kW -> kA conversion]
    Vtf = 8320                             # distr-level transformer rms voltage --- V [used in inelastic kW -> kA conv]
    Ntf   = Vtf/Vac                        # pole-top transformer turns ratio

    # Discretization parameters:
    #Ts = Rh*C/9              # s, sampling time in seconds
    Ts=180.0
    #ηP=a*(Ts/3600)*Vac #Vh
    ηP=a*(Ts/3600)*Vac/1e3 #kVh
    τP = exp(- Ts/(Rh*C))
    ρP = 1 - τP            # no units, ambient-to-temp param: 1/RC
    γP = Rh*Rw/(Ntf)*ρP*1000^2/1000    # kK/kW, ohmic losses-to-temp parameter

    ## MPC Paramters
    T1=12 #hours
    T2=14-T1
    # T1=6 #hours
    # T2=14-T1
    K1 = round(Int,T1*3600/Ts);            # Initial Prediction and Fixed Horizon (assume K1 instants = 12 hrs)
    K2 = round(Int,T2*3600/Ts);             # Additional time instants past control horizon
    K  = K1+K2;                        # Total horizon (8 PM to 10 AM)

    # PWL Parameters:
    S=15
    ItotalMax = 4  #kA
    deltaI = ItotalMax/S

    #system information
    t0=.370
    Tamb=.37*ones(K,1)
    iD_pred=0*ones(K,1)
    iD_actual=iD_pred

    Qmag=1
    Rmag=1
    Omag=100
    Rh=(Rmag*rand(1,H).+Rmag/1e3)
    Qh=(Qmag*rand(1,H).+Qmag/1e3)
    Oh=(Omag*rand(1,H).+Omag/1e3)

    #action happens immediately afte interval before ends
    #hub information
    arriveLast=round(Int,2*3600/Ts) #last arrival at 10PM
    departFirst=round(Int,10*3600/Ts) #first departure at 6AM
    K_arrive_pred=rand(1:arriveLast,Nh,H)
    K_depart_pred=rand(departFirst:K,Nh,H)
    # K_arrive_pred=zeros(Nh,H)
    # K_arrive_pred[1:Int(Nh/2),1]=1:Int(Nh/2)
    # K_depart_pred=hcat(1:Nh)
    # K_arrive_pred=hcat(1:Nh)
    # K_depart_pred=hcat(K-Nh:K)
    K_arrive_actual=K_arrive_pred
    K_depart_actual=K_depart_pred
    Sn_depart_min=round.(1 .- 0.20*rand(Nh,H),digits=4) #need 80-100%
    #Sn_depart_min=.8*ones(Nh,H)
    Sn_arrive_pred=round.(0.20*rand(Nh,H),digits=4) #arrive with 0-20%
    #Sn_arrive_pred=.8*ones(Nh,H)
    Sn_arrive_actual=Sn_arrive_pred
    #EVcap=b./3.6e6 #kWh
    EVcap=round.(b./3.6e6/1e3,digits=6) #MWh
    e0=zeros(H) # no vehicles have arrived yet
    #e0=[sum(Sn_arrive_actual[n,h]*EVcap[n,h] for n in findall(x->x==0,K_arrive_actual[:,h])) for h=1:H]

    #prepare predicted values for optimization
    eMax=zeros(K,H)
    uMax=zeros(K,H)
    eDepart_min=zeros(K,H)
    eArrive_pred=zeros(K,H)
    eArrive_actual=zeros(K,H)
    slackMax=zeros(K,H)
    for k =1:K
        for h=1:H # find a way to do this without another loop
            depart=[n for n=1:Nh if k==K_depart_pred[n,h]]
            if length(depart)!=0
                eDepart_min[k,h] = sum(Sn_depart_min[n,h]*EVcap[n,h] for n in depart)
                slackMax[k,h] = sum(EVcap[n,h]-Sn_depart_min[n,h]*EVcap[n,h] for n in depart)
            end

            # arrive=[n for n=1:Nh,h=1:H if k==hubS.K_arrive_pred[n,h]]
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

    hubS=scenarioHubStruct(Nh,H,Ts,K1,K2,K,S,ItotalMax,deltaI,Tmax,ηP,τP,ρP,γP,e0,t0,iD_pred,iD_actual,Tamb,Qh,Rh,Oh,
                            Sn_depart_min,Sn_arrive_actual,Sn_arrive_pred,K_arrive_pred,K_depart_pred,K_arrive_actual,
                            K_depart_actual,EVcap,eMax,uMax,eDepart_min,eArrive_pred,eArrive_actual,slackMax)

    if saveS==true
            save(path*"HubscenarioH$(H).jld2","hubS",hubS)
    end

    return hubS
end
