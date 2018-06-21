using Parameters

#structures
@with_kw struct convMetrics
    objVal::Array=zeros(maxIt,1)
    couplConst::Array=zeros(maxIt,1)
    lam::Array=zeros(maxIt,1)
    sn::Array=zeros(maxIt,1)
    un::Array=zeros(maxIt,1)

    lamIt::Array=zeros(maxIt,1)
    snIt::Array=zeros(maxIt,1)
    unIt::Array=zeros(maxIt,1)
end

@with_kw struct itLogPWL
    #model variables
    Xt::Array{Float64}=zeros((horzLen+1),maxIt) #rows are time
    Sn::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt)
    Un::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt)
    Z::Array{Float64}=zeros(S*(horzLen+1),maxIt)  #row are time,  columns are iteration

    #extra model variables
    Tactual::Array{Float64}=zeros((horzLen+1),maxIt) #rows are time
    uSum::Array{Float64}=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration
    zSum::Array{Float64}=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration
    couplConst::Array{Float64}=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration

    #dual vairables
    Lam::Array{Float64}=zeros((horzLen+1),maxIt) #(rows are time, columns are iteration)

    #auxillary variables
    Vu::Array{Float64}=zeros((N)*(horzLen+1),maxIt) #row are time,  columns are iteration
    Vz::Array{Float64}=zeros(S*(horzLen+1),maxIt)
    Vs::Array{Float64}=zeros((N)*(horzLen+1),maxIt) #row are time,  columns are iteration
    Vt::Array{Float64}=zeros((horzLen+1),maxIt) #row are time,  columns are iteration

    #Gradian Vectors
    Gu::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt) #row are time (N states for k=1, them N states for k=2),  columns are iteration
    Gs::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt) #row are time (N states for k=1, them N states for k=2),  columns are iteration
    Gz::Array{Float64}=zeros(S*(horzLen+1),maxIt) #row are time (N states for k=1, them N states for k=2),  columns are iteration

    #Jacobian C Vectors
    Cs::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt)  #row are time,  columns are iteration
    Cu::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt)  #row are time,  columns are iteration
    Cz::Array{Float64}=zeros(S*(horzLen+1),maxIt)  #row are time,  columns are iteration
    Ct::Array{Float64}=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration
end

@with_kw struct itLogNL
    #model variables
    Xt::Array{Float64}=zeros((horzLen+1),maxIt) #rows are time
    Sn::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt)
    Un::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt)
    Itotal::Array{Float64}=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration

    #extra model variables
    uSum::Array{Float64}=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration
    couplConst::Array{Float64}=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration

    #dual vairables
    Lam::Array{Float64}=zeros((horzLen+1),maxIt) #(rows are time, columns are iteration)

    #auxillary variables
    Vu::Array{Float64}=zeros((N)*(horzLen+1),maxIt) #row are time,  columns are iteration
    Vi::Array{Float64}=zeros((horzLen+1),maxIt)
    Vs::Array{Float64}=zeros((N)*(horzLen+1),maxIt) #row are time,  columns are iteration
    Vt::Array{Float64}=zeros((horzLen+1),maxIt) #row are time,  columns are iteration

    #Gradian Vectors
    Gu::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt) #row are time (N states for k=1, them N states for k=2),  columns are iteration
    Gs::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt) #row are time (N states for k=1, them N states for k=2),  columns are iteration
    Gi::Array{Float64}=zeros((horzLen+1),maxIt) #row are time (N states for k=1, them N states for k=2),  columns are iteration

    #Jacobian C Vectors
    Cs::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt)  #row are time,  columns are iteration
    Cu::SharedArray{Float64}=zeros(N*(horzLen+1),maxIt)  #row are time,  columns are iteration
    Ci::Array{Float64}=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration
    Ct::Array{Float64}=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration
end
#
# @with_kw struct centralSolution
#     xt::Array
#     sn::Array
#     un::Array
#     z::Array #for PWL
#     itotal::Array #for NL
#     objVal::Float64
#     lamTemp::Array
#     lamCoupl::Array
#     uSum::Array
#     zSum::Array
#     Tactual::Array
# end

@with_kw struct centralSolution
    xt::Array=zeros((horzLen+1),1)
    sn::Array=zeros(N*(horzLen+1),1)
    un::Array=zeros(N*(horzLen+1),1)
    z::Array=zeros(S*(horzLen+1),1) #for PWL
    itotal::Array=zeros((horzLen+1),1) #for NL
    objVal::Float64=0
    lamTemp::Array=zeros((horzLen+1),1)
    lamCoupl::Array=zeros((horzLen+1),1)
    uSum::Array=zeros(N*(horzLen+1),1)
    zSum::Array=zeros(S*(horzLen+1),1)
    Tactual::Array=zeros((horzLen+1),1)
end
