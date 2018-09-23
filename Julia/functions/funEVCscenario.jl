#generate simulation scenario for Centralized EV Charging Problem
#Micah Botkin-Levy
#4/13/18
#from Mads Almassakhi code

#structures
struct scenarioStruct
    N::Int

    #horizon
    K1::Int
    K2::Int
    K::Int

    #PWL
    S::Int
    ItotalMax::Int
    deltaI::Float64

    #limits
    Tmax::Float64
    imin::Array{Float64,2} #switch these to 1 dim array/vectors
    imax::Array{Float64,2}

    #Discretization Paramters
    ηP::Array{Float64,2}
    τP::Float64
    ρP::Float64
    γP::Float64

    #initial conditions
    s0::Array{Float64,2}
    t0::Int

    #desired conditions
    Snmin::Array{Float64,2}
    Kn::Array{Int,2}

    #disturbances
    w::Array{Float64,2}

    #User def penalty matrix
    Qsi::Array{Float64,2}
    Ri::Array{Float64,2}
end

#functions
function setupScenario(N;Tmax=393,Dload_amplitude=0,saveS=false)

    #model parameters
    a   = rand(N,1)*.1 + 0.8               # efficiency of Li-ion batts is ~80-90%
    b   = (6*rand(N,1)+12)*3.6e6           # battery capacity (12-18 kWh = 43.3-64.8 MJ)
    m   = 2000                             # transformer mass in kg
    C   = 450*m                            # heat cap. thermal mass J/K ----- spec. heat cap. of C = {carbon steel, iron, veg. oil} = {490, 450, 1670} J/(kg*K)
    Rh   = 1070e-4/(35*5*(m/7870)^(2/3))   # heat outflow resistance K/W : R = 0.1073 (K*m^2/W)/(A_s), rule of thumb calculation
    Rw  = 1                                # coil winding resistance --- ohms:
    Vac = 240                              # PEV battery rms voltage --- V [used in PEV kW -> kA conversion]
    Vtf = 8320                             # distr-level transformer rms voltage --- V [used in inelastic kW -> kA conv]
    Ntf   = Vtf/Vac                        # pole-top transformer turns ratio

    # Discretization parameters:
    Ts = Rh*C/9              # s, sampling time in seconds
    ηP = Ts*Vac*a./b*1000  # 1/kA, normalized battery sizes (0-1)
    #ηP = Ts*Vac*a./b  # 1/A, normalized battery sizes (0-1)
    τP = 1 - Ts/(Rh*C)      # no units, temp time constant: 1 - 1/RC
    ρP = 1 - τP            # no units, ambient-to-temp param: 1/RC
    γP = Ts*Rw/(C*Ntf)*1000^2    # K/kW, ohmic losses-to-temp parameter
    #γP = Ts*Rw/(C*Ntf)    # K/W, ohmic losses-to-temp parameter

    # PWL Parameters:
    #S = 3;
    S=10
    #ItotalMax = 20;        % CAUTION  ---> Imax gives upper limit on total current input on Transfomer and if picked too low will cause infeasible.
    ItotalMax = 4  #kA
    #ItotalMax = 4000  #A
    deltaI = ItotalMax/S


    ## MPC Paramters
    T1=12
    T2=2
    K1 = round(Int,T1*3600/Ts);            # Initial Prediction and Fixed Horizon (assume K1 instants = 12 hrs)
    K2 = round(Int,T2*3600/Ts);             # Additional time instants past control horizon
    K  = K1+K2;                        # Total horizon (8 PM to 10 AM)Qs  = 10;              % Stage and terminal penalty on charge difference with respect to 1 (states s)
    #K=331
    #K=199

    # Constraint parameters:
    #Tmax = 393                             # Short-term over-loading --> 120 C = 393 Kelvin
    imin = zeros(N,1)                      # A, q_min < 0 if V2G is allowed
    imax = (10 + 16*rand(N,1))/1000             # kA, charging with 10-24 A
    #imax = (10 + 16*rand(N,1))             # A, charging with 10-24 A

    # Initial conditions:
    s0 = 0.2*rand(N,1)      # initial states of charge (0 - 0.20)
    t0 = 370                 # initial temp (~65 K below Tmax) 368K

    #desired states
    SOCmin = 1 - 0.20*rand(N,1)            # Required min final states of charge (~0.80-1)
    FullChargeTime_relative = .25*rand(N,1)+.75
    FullChargeTime = convert(Array{Int,2},round.(K1*FullChargeTime_relative))
    Snmin=SOCmin
    Kn=FullChargeTime

    # Disturbances
    #Dload_amplitude = 2;  # base-demand factor
    #Dload_amplitude = 85 #kWatts?
    #Dload_amplitude = 75000 #Watts?
    #Dload_amplitude = 0
    Tamb_amplitude  = 370   # assume hot night in summer (30 C) 363K

    # Disturbance scenario:
    #FullinelasticDemand = [normpdf(0,linspace(0,8,round((K1-1)/2)),3) normpdf(0,linspace(-8,0,round(K1/2)),3)]; # let demand per household be peaking at 8PM and 8 PM with nadir inbetween
    #FullinelasticDemand = 100*(200*(FullinelasticDemand-min(FullinelasticDemand))/range(FullinelasticDemand) + 600)/1000; # total non-EV demand (in kW) = N/PEVpenetration*inelasticDemandperHouse
    #FullinelasticDemand = [FullinelasticDemand'; FullinelasticDemand(end)*ones(K2+1,1)];
    #inelasticDemand = [normpdf(0,linspace(0,8,round(K/2)),3) normpdf(0,linspace(-8,0,round(K/2)),3)]; # let demand per household be peaking at 8PM and 8 PM with nadir inbetween

    #dist=[linspace(0,8,round(K/2));linspace(-8,0,round(K/2))]
    dist = [linspace(0,8,round(K1/2));linspace(-8,0,K1-round(K1/2))]
    d = Normal(0,3)
    inelasticDemand = pdf.(d,dist)
    FullinelasticDemand = 100*(200*(inelasticDemand-minimum(inelasticDemand))/(maximum(inelasticDemand)-minimum(inelasticDemand)) + 600)/1000; # total non-EV demand (in kW) = N/PEVpenetration*inelasticDemandperHouse
    FullinelasticDemand = [FullinelasticDemand; FullinelasticDemand[length(FullinelasticDemand)]*ones(K2+1,1)]
    FullDload   = Dload_amplitude*FullinelasticDemand;    # peaks during mid-day
    iD = FullDload/Vtf;                                     #background demand current
    FullTamb    = Tamb_amplitude*ones(K+1,1);             #normpdf(0,linspace(-10,10,max(K,kmax)),3)';   # exogenous peaks during mid-day          % OVER-NIGHT CHARGING: TIMES -1?
    ndisturbs = 2;
    w = zeros((K+1)*ndisturbs,1);
    for i=1:K+1
        w[(i-1)*ndisturbs+1:i*ndisturbs,1]  = [iD[i,1]; FullTamb[i,1]];
    end

    # penalty matrix new (need to fix for k>Ki)
    Ru   = 0.1*1000^2;              # Stage and terminal penalty on local power flow (inputs u)
    #Ru   = 0.1;              # Stage and terminal penalty on local power flow (inputs u)
    #RKi   = 10;            # Stage and terminal penalty on local power flow (inputs q), for k >= Ki
    Qs  = 10;               # Stage and terminal penalty on charge difference with respect to 1 (states s)
    QT  = 0;                # PENALTY ON TEMPERATURE DEVIATION (W.R.T 0)
    #QsKi  = 1;             # Stage and terminal penalty on charge difference with respect to 1 (states s), for k >= Ki
    Ri=Ru*(5*rand(N,1)+.1);
    Qsi=[Qs*(10*rand(N,1)+.01);QT];

    #move this into struct???
    @assert all(ηP.*K.*FullChargeTime_relative.*imax+s0 .>= SOCmin) "Some PEVs may not be able to meet SOC min level by desired time!"


    evScenario=scenarioStruct(N,K1,K2,K,S,ItotalMax,deltaI,Tmax,imin,imax,
                            ηP,τP,ρP,γP,s0,t0,Snmin,Kn,w,Qsi,Ri)

    if saveS==true
        JLD.save("EVCscenarioN$(N).jld","evScenario",evScenario)
    end

    return evScenario
end