#comparison plots


#compare convergence rates
legendLab=["Dual Ascent", "Fast Ascent", "ADMM"]
legendCol=["blue","red", "green"]

#f-fStar
fgapPlotComp=plot(layer(x=1:convIt,y=fConvDual_dualAscent[1:convIt,1],Geom.line,Theme(default_color=colorant"blue")),
            layer(x=1:convIt,y=fConvDual_fastAscent[1:convIt,1],Geom.line,Theme(default_color=colorant"red")),
            layer(x=1:convIt,y=fConvADMM[1:convIt,1],Geom.line,Theme(default_color=colorant"green")),
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Objective Function Gap",orientation=:vertical),Scale.y_log10,
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt),
            Guide.manual_color_key("", legendLab, legendCol))

#norm(u-uStar)
unPlotComp=plot(layer(x=1:convIt,y=unConvDual_dualAscent[1:convIt,1],Geom.line,Theme(default_color=colorant"blue")),
            layer(x=1:convIt,y=unConvDual_fastAscent[1:convIt,1],Geom.line,Theme(default_color=colorant"red")),
            layer(x=1:convIt,y=unConvADMM[1:convIt,1],Geom.line,Theme(default_color=colorant"green")),
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm EV Current Gap",orientation=:vertical),Scale.y_log10,
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt),
            Guide.manual_color_key("", legendLab, legendCol))


#Ax-b=0
constPlotComp=plot(layer(x=1:convIt,y=constConvDual[1:convIt,1],Geom.line,Theme(default_color=colorant"blue")),
            layer(x=1:convIt,y=constConvADMM[1:convIt,1],Geom.line,Theme(default_color=colorant"green")),
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Ax Gap",orientation=:vertical),Scale.y_log10,
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt),
            Guide.manual_color_key("", ["Dual Ascent", "ADMM"], ["blue", "green"]))

#lambda-lamstar
lamPlotComp=plot(layer(x=1:convIt,y=ConvDual_dualAscent[1:convIt,1],Geom.line,Theme(default_color=colorant"blue")),
            layer(x=1:convIt,y=ConvDual_fastAscent[1:convIt,1],Geom.line,Theme(default_color=colorant"red")),
            layer(x=1:convIt,y=ConvADMM[1:convIt,1],Geom.line,Theme(default_color=colorant"green")),
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap",orientation=:vertical),Scale.y_log10,
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt),
            Guide.manual_color_key("", legendLab, legendCol))

#lambda-lam[p-1]
itPlotComp=plot(layer(x=1:convIt,y=itConvDual_dualAscent[1:convIt,1],Geom.line,Theme(default_color=colorant"blue")),
            layer(x=1:convIt,y=itConvDual_fastAscent[1:convIt,1],Geom.line,Theme(default_color=colorant"red")),
            layer(x=1:convIt,y=itConvADMM[1:convIt,1],Geom.line,Theme(default_color=colorant"green")),
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Iterate Gap",orientation=:vertical),Scale.y_log10,
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt),
            Guide.manual_color_key("", legendLab, legendCol))


path="C:\\Users\\micah\\Documents\\uvm\\Research\\EVC code\\Julia\\convergence plots\\"

if drawFig==1 draw(PNG(path*"itconv_N$(N).png", 36inch, 24inch), itPlotComp) end


if drawFig==1 draw(PNG(path*"conv_N$(N).png", 36inch, 24inch), hstack(vstack(fgapPlotComp,unPlotComp),vstack(constPlotComp,lamPlotComp))) end
