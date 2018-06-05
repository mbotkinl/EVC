#comparison plots


#compare convergence rates

#f-fStar
fgapPlotComp=plot(layer(x=1:convIt,y=fConvDual[1:convIt,1],Geom.line,Theme(default_color=colorant"blue")),
            layer(x=1:convIt,y=fConvADMM[1:convIt,1],Geom.line,Theme(default_color=colorant"green")),
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Objective Function Gap",orientation=:vertical),Scale.y_log10,
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt),
            Guide.manual_color_key("", ["Dual Ascent", "ADMM"], ["blue", "green"]))

#norm(u-uStar)
unPlotComp=plot(layer(x=1:convIt,y=unConvDual[1:convIt,1],Geom.line,Theme(default_color=colorant"blue")),
            layer(x=1:convIt,y=unConvADMM[1:convIt,1],Geom.line,Theme(default_color=colorant"green")),
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm EV Current Gap",orientation=:vertical),Scale.y_log10,
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt),
            Guide.manual_color_key("", ["Dual Ascent", "ADMM"], ["blue", "green"]))


#Ax-b=0
constPlotComp=plot(layer(x=1:convIt,y=constConvDual[1:convIt,1],Geom.line,Theme(default_color=colorant"blue")),
            layer(x=1:convIt,y=constConvADMM[1:convIt,1],Geom.line,Theme(default_color=colorant"green")),
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Ax Gap",orientation=:vertical),Scale.y_log10,
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt),
            Guide.manual_color_key("", ["Dual Ascent", "ADMM"], ["blue", "green"]))

#lambda-lamstar
lamPlotComp=plot(layer(x=1:convIt,y=ConvDual[1:convIt,1],Geom.line,Theme(default_color=colorant"blue")),
            layer(x=1:convIt,y=ConvADMM[1:convIt,1],Geom.line,Theme(default_color=colorant"green")),
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap",orientation=:vertical),Scale.y_log10,
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt),
            Guide.manual_color_key("", ["Dual Ascent", "ADMM"], ["blue", "green"]))

path="C:\\Users\\micah\\Documents\\uvm\\Research\\EVC code\\Julia\\convergence plots\\"
if drawFig==1 draw(PNG(path*"conv_$(N).png", 36inch, 12inch), hstack(vstack(fgapPlotComp,unPlotComp),vstack(constPlotComp,lamPlotComp))) end
