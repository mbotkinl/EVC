

#scenarios
struct scenarioStruct
    N::Int
    Ts::Float64

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
    t0::Float64

    #desired conditions
    Snmin::Array{Float64,2}
    Kn::Array{Int,2}

    #disturbances
    iD_pred::Array{Float64,2}
    iD_actual::Array{Float64,2}
    Tamb::Array{Float64,2}

    #User def penalty matrix
    Qsi::Array{Float64,2}
    Ri::Array{Float64,2}
    β::Array{Float64,2}
end

struct scenarioHubStruct
    Nh::Int64
    H::Int64
    Ts::Float64

    #horizon
    K1::Int64
    K2::Int64
    K::Int64

    #PWL
    S::Int64
    ItotalMax::Float64
    deltaI::Float64

    #limits
    Tmax::Float64

    #Discretization Paramters
    ηP::Array{Float64,2}
    τP::Float64
    ρP::Float64
    γP::Float64

    #initial conditions
    e0::Array{Float64,1}
    t0::Float64

    #disturbances
    iD_pred::Array{Float64,2}
    iD_actual::Array{Float64,2}
    Tamb::Array{Float64,2}

    #User def penalty matrix
    Qh::Array{Float64,2}
    Rh::Array{Float64,2}
    Oh::Array{Float64,2}

    #hub conditions
    Sn_depart_min::Array{Float64,2}
    Sn_arrive_actual::Array{Float64,2}
    Sn_arrive_pred::Array{Float64,2}
    K_arrive_pred::Array{Int64,2}
    K_depart_pred::Array{Int64,2}
    K_arrive_actual::Array{Int64,2}
    K_depart_actual::Array{Int64,2}
    EVcap::Array{Float64,2}

    eMax::Array{Float64,2}
    uMax::Array{Float64,2}
    eDepart_min::Array{Float64,2}
    eArrive_pred::Array{Float64,2}
    eArrive_actual::Array{Float64,2}
    slackMax::Array{Float64,2}
end

#debug info
@with_kw struct convMetricsStruct
    obj::Array=zeros(maxIt,1)
    couplConst::Array=zeros(maxIt,1)
    lam::Array=zeros(maxIt,1)
    sn::Array=zeros(maxIt,1)
    un::Array=zeros(maxIt,1)

    lamIt::Array=zeros(maxIt,1)
    snIt::Array=zeros(maxIt,1)
    unIt::Array=zeros(maxIt,1)
end

@with_kw struct itLogPWL
    horzLen::Int
    N::Int
    S::Int
    # objective value
    objVal::Array{Float64}=zeros(1,maxIt+1) #columns are iteration

    #model variables
    Tpwl::Array{Float64}=zeros((horzLen+1),maxIt+1) #rows are time
    Sn::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt+1)
    Un::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt+1)
    Z::Array{Float64}=zeros(S*(horzLen+1),maxIt+1)  #row are time,  columns are iteration
    slackSn::Array=zeros(N,1)

    #extra model variables
    Tactual::Array{Float64}=zeros((horzLen+1),maxIt+1) #rows are time
    Itotal::Array{Float64}=zeros((horzLen+1),maxIt+1) #rows are time
    uSum::Array{Float64}=zeros((horzLen+1),maxIt+1)  #row are time,  columns are iteration
    zSum::Array{Float64}=zeros((horzLen+1),maxIt+1)  #row are time,  columns are iteration
    couplConst::Array{Float64}=zeros((horzLen+1),maxIt+1)  #row are time,  columns are iteration

    #dual vairables
    Lam::Array{Float64}=zeros((horzLen+1),maxIt+1) #(rows are time, columns are iteration)

    #auxillary variables
    Vu::Array{Float64}=zeros((N)*(horzLen+1),maxIt+1) #row are time,  columns are iteration
    Vz::Array{Float64}=zeros(S*(horzLen+1),maxIt+1)
    Vs::Array{Float64}=zeros((N)*(horzLen+1),maxIt+1) #row are time,  columns are iteration
    Vt::Array{Float64}=zeros((horzLen+1),maxIt+1) #row are time,  columns are iteration

    #Gradian Vectors
    Gu::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt+1) #row are time (N states for k=1, them N states for k=2),  columns are iteration
    Gs::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt+1) #row are time (N states for k=1, them N states for k=2),  columns are iteration
    Gz::Array{Float64}=zeros(S*(horzLen+1),maxIt+1) #row are time (N states for k=1, them N states for k=2),  columns are iteration
    Gt::Array{Float64}=zeros((horzLen+1),maxIt+1) #row are time (N states for k=1, them N states for k=2),  columns are iteration

    #Jacobian C Vectors
    Csu::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt+1)  #row are time,  columns are iteration
    Cuu::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt+1)  #row are time,  columns are iteration
    Czu::Array{Float64}=zeros(S*(horzLen+1),maxIt+1)  #row are time,  columns are iteration
    Ctu::Array{Float64}=zeros((horzLen+1),maxIt+1)  #row are time,  columns are iteration
    Csl::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt+1)  #row are time,  columns are iteration
    Cul::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt+1)  #row are time,  columns are iteration
    Czl::Array{Float64}=zeros(S*(horzLen+1),maxIt+1)  #row are time,  columns are iteration
    Ctl::Array{Float64}=zeros((horzLen+1),maxIt+1)  #row are time,  columns are iteration

    #update dynamic (either rho or alpha)
    itUpdate::Array{Float64}=zeros(1,maxIt+1) #columns are iteration
end

@with_kw struct itLogNL
    # objective value
    objVal::Array{Float64}=zeros(1,maxIt+1) #columns are iteration

    #model variables
    Xt::Array{Float64}=zeros((horzLen+1),maxIt+1) #rows are time
    Sn::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt+1)
    Un::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt+1)
    Itotal::Array{Float64}=zeros((horzLen+1),maxIt+1)  #row are time,  columns are iteration
    slackSn::Array=zeros(N,1)

    #extra model variables
    Tactual::Array{Float64}=zeros((horzLen+1),maxIt+1) #rows are time
    uSum::Array{Float64}=zeros((horzLen+1),maxIt+1)  #row are time,  columns are iteration
    couplConst::Array{Float64}=zeros((horzLen+1),maxIt+1)  #row are time,  columns are iteration

    #dual vairables
    Lam::Array{Float64}=zeros((horzLen+1),maxIt+1) #(rows are time, columns are iteration)

    #auxillary variables
    Vu::Array{Float64}=zeros((N)*(horzLen+1),maxIt+1) #row are time,  columns are iteration
    Vi::Array{Float64}=zeros((horzLen+1),maxIt+1)
    Vs::Array{Float64}=zeros((N)*(horzLen+1),maxIt+1) #row are time,  columns are iteration
    Vt::Array{Float64}=zeros((horzLen+1),maxIt+1) #row are time,  columns are iteration

    #Gradian Vectors
    Gu::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt+1) #row are time (N states for k=1, them N states for k=2),  columns are iteration
    Gs::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt+1) #row are time (N states for k=1, them N states for k=2),  columns are iteration
    Gi::Array{Float64}=zeros((horzLen+1),maxIt+1) #row are time (N states for k=1, them N states for k=2),  columns are iteration
    Gt::Array{Float64}=zeros((horzLen+1),maxIt+1) #row are time (N states for k=1, them N states for k=2),  columns are iteration

    #Jacobian C Vectors
    Csu::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt+1)  #row are time,  columns are iteration
    Cuu::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt+1)  #row are time,  columns are iteration
    Ciu::Array{Float64}=zeros((horzLen+1),maxIt+1)  #row are time,  columns are iteration
    Ctu::Array{Float64}=zeros((horzLen+1),maxIt+1)  #row are time,  columns are iteration
    Csl::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt+1)  #row are time,  columns are iteration
    Cul::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt+1)  #row are time,  columns are iteration
    Cil::Array{Float64}=zeros((horzLen+1),maxIt+1)  #row are time,  columns are iteration
    Ctl::Array{Float64}=zeros((horzLen+1),maxIt+1)  #row are time,  columns are iteration

    #update dynamic (either rho or alpha)
    itUpdate::Array{Float64}=zeros(1,maxIt+1) #columns are iteration
end

#solutions
@with_kw struct solutionStruct
    K::Int
    N::Int
    S::Int
    Tpwl::Array=zeros(K,1)
    Sn::Array=zeros(K,N)
    Un::Array=zeros(K,N)
    Z::Array=zeros(K,S) #for PWL only
    Itotal::Array=zeros(K,1)
    objVal::Array=zeros(1,1)
    lamTemp::Array=zeros(K,1)
    lamCoupl::Array=zeros(K,1)
    uSum::Array=zeros(K,1)
    zSum::Array=zeros(K,1)
    Tactual::Array=zeros(K,1)
    slackSn::Array=zeros(K,1)
    convIt::Array=zeros(K,1)
end

@with_kw struct hubSolutionStruct
    K::Int
    H::Int
    E::Array=zeros(K,H) #row are time, column are hub
    U::Array=zeros(K,H) #row are time, column are hub
    Itotal::Array=zeros(K,1) #row are time
    Tpwl::Array=zeros(K,1) #row are time
    Tactual::Array=zeros(K,1) #row are time
    uSum::Array=zeros(K,1) #row are time
    zSum::Array=zeros(K,1) #row are time
    Lam::Array=zeros(K,1) #row are time
    convIt::Array=zeros(K,1) #row are time
    objVal::Array=zeros(K,1)

    E_depart::Array=zeros(K,H)
    E_arrive::Array=zeros(K,H)
end
#
# @with_kw struct hubItAux
#     horzLen::Int
#     H::Int
#     S::Int
#     K::Int
#     T0::Array=zeros(K,1)  #row are time, column are hub
#     E0::Array=zeros(K,H)  #row are time, column are hub
#     Ve::Array=zeros(horzLen+1,H,maxIt+1) #row are time, column are hub, stack is iteration
#     Vu::Array=zeros(horzLen+1,H,maxIt+1) #row are time, column are hub, stack is iteration
#     Vz::Array=zeros(horzLen+1,S,maxIt+1) #row are time, column are segment, stack is iteration
#     Vt::Array=zeros(horzLen+1,1,maxIt+1) #row are time, stack is iteration
#     VΔ::Array=zeros(horzLen+1,H,maxIt+1) #row are time, stack is iteration
#     Lam::Array=zeros(horzLen+1,1,maxIt+1) #row are time, stack is iteration
# end

@with_kw struct hubItLogPWL
    horzLen::Int
    H::Int
    S::Int
    E::Array{Float64,3}=zeros(horzLen+1,H,maxIt) #row are time, column are hub, stack is iteration
    U::Array{Float64,3}=zeros(horzLen+1,H,maxIt) #row are time, column are hub, stack is iteration
    Z::Array{Float64,3}=zeros(horzLen+1,S,maxIt) #row are time, column are segment, stack is iteration
    Tpwl::Array{Float64,3}=zeros(horzLen+1,1,maxIt) #row are time, stack is iteration
    Lam::Array{Float64,3}=zeros(horzLen+1,1,maxIt) #row are time, stack is iteration
    couplConst::Array{Float64,3}=zeros(horzLen+1,1,maxIt) #row are time, stack is iteration
    Itotal::Array{Float64,3}=zeros(horzLen+1,1,maxIt) #row are time, stack is iteration
    uSum::Array{Float64,3}=zeros((horzLen+1),1,maxIt+1) #row are time, stack is iteration
    zSum::Array{Float64,3}=zeros((horzLen+1),1,maxIt+1) #row are time, stack is iteration
    Tactual::Array{Float64,3}=zeros(horzLen+1,1,maxIt) #row are time, stack is iteration

    #auxillary variables
    Vu::Array{Float64,3}=zeros((horzLen+1),H,maxIt+1)#row are time, columns are hub,  stacks are iteration
    Vz::Array{Float64,3}=zeros((horzLen+1),S,maxIt+1)
    Ve::Array{Float64,3}=zeros((horzLen+1),H,maxIt+1) #row are time, columns are hub,  stacks are iteration
    Vd::Array{Float64,3}=zeros((horzLen+1),H,maxIt+1) #row are time, columns are hub,  stacks are iteration
    Vt::Array{Float64,3}=zeros((horzLen+1),1,maxIt+1) #row are time,  columns are iteration

    #Gradian Vectors
    Gu::SharedArray{Float64,3}=zeros((horzLen+1),H,maxIt+1) #row are time, columns are hub,  stacks are iteration
    Ge::SharedArray{Float64,3}=zeros((horzLen+1),H,maxIt+1) #row are time, columns are hub,  stacks are iteration
    GeΔ::SharedArray{Float64,3}=zeros((horzLen+1),H,maxIt+1) #row are time, columns are hub,  stacks are iteration
    Gz::Array{Float64,3}=zeros((horzLen+1),S,maxIt+1) #row are time, stacks are iteration
    Gt::Array{Float64,3}=zeros((horzLen+1),1,maxIt+1) #row are time,   stacks are iteration

    #Jacobian C Vectors
    Ceu::SharedArray{Float64,3}=zeros((horzLen+1),H,maxIt+1)   #row are time, columns are hub,  stacks are iteration
    Cuu::SharedArray{Float64,3}=zeros((horzLen+1),H,maxIt+1)   #row are time, columns are hub,  stacks are iteration
    Cdu::SharedArray{Float64,3}=zeros((horzLen+1),H,maxIt+1)   #row are time, columns are hub,  stacks are iteration
    Czu::Array{Float64,3}=zeros((horzLen+1),S,maxIt+1)  #row are time,   stacks are iteration
    Ctu::Array{Float64,3}=zeros((horzLen+1),1,maxIt+1)  #row are time,   stacks are iteration

    Cel::SharedArray{Float64,3}=zeros((horzLen+1),H,maxIt+1)   #row are time, columns are hub,  stacks are iteration
    Cul::SharedArray{Float64,3}=zeros((horzLen+1),H,maxIt+1)   #row are time, columns are hub,  stacks are iteration
    Cdl::SharedArray{Float64,3}=zeros((horzLen+1),H,maxIt+1)   #row are time, columns are hub,  stacks are iteration
    Czl::Array{Float64,3}=zeros((horzLen+1),S,maxIt+1)  #row are time,   stacks are iteration
    Ctl::Array{Float64,3}=zeros((horzLen+1),1,maxIt+1)  #row are time,   stacks are iteration

    itUpdate::Array{Float64,3}=zeros(1,1,maxIt+1)  #row are time,   stacks are iteration
end
