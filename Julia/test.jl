using Plots; gr()
n = 10; x = range(0,step=1,length=n); y=cumsum(rand(n));
plot(x,y,line_z=1:n, w=5)


# testing gradient plot
N=40
using Plots; pyplot()
n = 10; x = range(0,step=1,length=n); y=rand(n,N);
plot(x,y,line_z=1:n, w=5)

colors=[:red :green :blue]

g=cgrad(:blues)
colors=[g[1] g[10]  g[30]]
colors=[g[z] for z=range(1,stop=length(g.values),length=N)]
plot(x,y,seriescolor=colors)
plot(x,y,palette=:blues)

levels=range(0,stop=1000,length=N)
plot(x,y,seriescolor=ColorGradient([:red, :yellow]))

scatter(x,y,zcolor=1:N,c=cgrad(:inferno))
plot(x,y,line_z=1:N,seriescolor=cgrad(:inferno))
plot(x,y,c=ColorGradient(:blues))

A = [i for i=1:100, j=1:100]
heatmap(A, c=ColorGradient([:red,:yellow,:blue]))
heatmap(A, c=cgrad(:inferno))

# plot(1:horzLen+1,dLog.uSum[:,1:4], linecolor=cgrad(:inferno))
# plot(1:horzLen+1,dLog.uSum[:,1:10],	color=cgrad(:grays))
# plot(1:horzLen+1,dLog.uSum[:,1:10],	zcolor=1:10)
#
# plot(1:horzLen+1,dLog.uSum[:,1:10], palette=:blues)


# C(g::ColorGradient) = RGB[g[z] for z=range(1,step=1,length=10)]
# g = :inferno
# cgrad(g) |> C
