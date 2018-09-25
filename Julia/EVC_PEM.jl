using Distributions
#pull out a few key variables
N=evS.N
S=evS.S
horzLen=evS.K1

#check indexes***
sn=zeros(horzLen+1,N)
un=zeros(horzLen+1,N)
R=zeros(horzLen,N)

sn[1:N]=evS.s0
for k =1:horzLen
    for n =1:N
        if sn[k,n]<evS.Snmin[n]
            R[k,n]=(evS.Snmin[n]-sn[k,n])/(evS.ηP[n]*evS.imax[n]*(evS.Kn[n]-k))
        else
            R[k,n]=0 # not sure if this is what we want***
        end
        t=rand()
        un[k+1,n]= if t>(1-R[k,n]) evS.imax[n] else  0  end #R[n,k]==0 #R[n,k]>=1
        sn[k+1,n]=sn[k,n]+evS.ηP[n]*un[k+1,n]
    end
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

pR=plot(R,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("R"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=24pt,
		minor_label_font_size=20pt,line_width=3pt,key_label_font_size=20pt))
