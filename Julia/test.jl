using Plots; gr()
n = 10; x = range(0,step=1,length=n); y=cumsum(rand(n));
plot(x,y,line_z=1:n, w=5)
