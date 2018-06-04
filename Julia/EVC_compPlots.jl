#comparison plots


#compare convergence rates

#f-fStar

convItPlotComp=plot(layer(x=1:convIt,y=itConvDual[1:convIt,1],Geom.line,Scale.y_log10,Theme(default_color=colorant"blue")),
            layer(x=1:convIt,y=itConvADMM[1:convIt,1],Geom.line,Scale.y_log10,Theme(default_color=colorant"green"))
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt),
            Guide.manual_color_key("", ["Dual Ascent", "ADMM"], ["blue", "green"]))
#norm(u-uStar)
