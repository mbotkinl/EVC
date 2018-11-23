using Parameters
using SharedArrays


#structures
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
    iD::Array{Float64,2}
    iDnoise::Array{Float64,2}
    Tamb::Array{Float64,2}

    #User def penalty matrix
    Qsi::Array{Float64,2}
    Ri::Array{Float64,2}
    β::Array{Float64,2}

end


#structures
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
    # objective value
    objVal::Array{Float64}=zeros(1,maxIt+1) #columns are iteration

    #model variables
    Xt::Array{Float64}=zeros((horzLen+1),maxIt+1) #rows are time
    Sn::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt+1)
    Un::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt+1)
    Z::Array{Float64}=zeros(S*(horzLen+1),maxIt+1)  #row are time,  columns are iteration
    slackSn::Array=zeros(N,1)

    #extra model variables
    Tactual::Array{Float64}=zeros((horzLen+1),maxIt+1) #rows are time
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

@with_kw struct centralSolutionStruct

    Xt::Array=zeros((horzLen+1),1)
    Sn::Array=zeros(N*(horzLen+1),1)
    Un::Array=zeros(N*(horzLen+1),1)
    z::Array=zeros(S*(horzLen+1),1) #for PWL
    Itotal::Array=zeros((horzLen+1),1) #for NL
    objVal::Float64=0
    lamTemp::Array=zeros((horzLen+1),1)
    lamCoupl::Array=zeros((horzLen+1),1)
    uSum::Array=zeros(N*(horzLen+1),1)
    zSum::Array=zeros(S*(horzLen+1),1)
    Tactual::Array=zeros((horzLen+1),1)
    slackSn::Array=zeros(N,1)
end
