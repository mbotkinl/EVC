using Distributions
#pull out a few key variables
N=evS.N
S=evS.S
horzLen=evS.K1

#check indexes***
sn=zeros(horzLen+1,N)
un=zeros(horzLen+1,N)
R=zeros(horzLen,N)
T=zeros(horzLen+1,1)

sn[1:N]=evS.s0
T[1]=evS.t0

for k =1:horzLen
    for n =1:N
        if sn[k,n]<evS.Snmin[n]
            R[k,n]=(evS.Snmin[n]-sn[k,n])/(evS.ηP[n]*evS.imax[n]*(evS.Kn[n]-k))
        else
            R[k,n]=0 # not sure if this is what we want***
        end
	end

	Iest = sum(if R[k,n]>0 evS.imax[n] else 0 end for n=1:N)
	Test = evS.τP*T[k]+evS.γP*Iest^2+evS.ρP*evS.w[k*2+2,1]

	for n=1:N
		# use probability
        t=rand()
        # un[k+1,n]= if t>(1-R[k,n]) evS.imax[n] else  0  end

		# decide only on temp
		un[k+1,n]= if R[k,n]>=1 evS.imax[n] elseif (Test<evS.Tmax) & (R[k,n]>0) evS.imax[n] else  0  end 

		#decided by Temp and probability
		un[k+1,n]= if R[k,n]>=1 evS.imax[n] elseif (Test<evS.Tmax) & (t>(1-R[k,n])) evS.imax[n] else  0  end

        sn[k+1,n]=sn[k,n]+evS.ηP[n]*un[k+1,n]
	end
	T[k+1] = evS.τP*T[k]+evS.γP*sum(un[k+1,n] for n=1:n)^2+evS.ρP*evS.w[k*2+2,1]
end


p1=plot(sn,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV SOC"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1,ymax=1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=24pt,line_width=3pt,
		minor_label_font_size=20pt,key_label_font_size=20pt))

p2=plot(un,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV Current (kA)"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=24pt,
		minor_label_font_size=20pt,line_width=3pt,key_label_font_size=20pt))

p3=plot(x=1:horzLen+1,y=T,Geom.line,Theme(default_color=colorant"green"),
		yintercept=[evS.Tmax],Geom.hline(color=["red"],style=:dot),
		Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",key_position = :top,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))

pR=plot(R,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("R"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=24pt,
		minor_label_font_size=20pt,line_width=3pt,key_label_font_size=20pt))
