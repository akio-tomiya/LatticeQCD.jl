module Realtimeplot
using Plots 

# prep the plots
plt1 = histogram([0],label=nothing) #plot1
plt2 = plot([], [],label=nothing) #plot2
ylabel!("Plaquette")
xlabel!("MC time")
plt3 = histogram([0],label=nothing) #plot3
plt4 = plot([], [],label=nothing) #plot4
ylabel!("Polyakov loop")
xlabel!("MC time")
plot(plt1, plt2, plt3, plt4, layout = 4)

function mock() #mock data
    plaq = 0.6+0.1*randn() 
    if rand()<0.5
        poly = -1.5+randn()/2
    else
        poly = 1.5+randn()/2
    end
    sleep(0.01)
    return plaq,poly
end

hist_plaq=[]
hist_poly=[]

function plot_refresh!(plt1,plt2,plt3,plt4,hist_plaq,hist_poly, plaq, poly,itrj)
    bins=round(Int,log(itrj)*5+1)
    append!(hist_plaq, plaq)
    append!(hist_poly, poly)
    #
    plt1 = histogram(hist_plaq,bins=bins,label=nothing) #plot1
    xlabel!("Plaquette")
    plt3 = histogram(hist_poly,bins=bins,label=nothing) #plot3
    xlabel!("Polyakov loop")
    #
    push!(plt2, 1, itrj, plaq)
    push!(plt4, 1, itrj, poly)
    #
    plot(plt1, plt2, plt3, plt4, layout = 4)
    gui()
end

Ntrj=1000
for itrj in 1:Ntrj
    plaq,poly=mock()
    #
    plot_refresh!(plt1,plt2,plt3,plt4,hist_plaq,hist_poly, plaq, poly,itrj)
    #
end
#=
plt = plot([0,0.1], Any[rand(2),sin])
for x in 0.2:0.1:π
    push!(plt, 1, x, rand())
    push!(plt, 2, x, sin(x))
    gui(); sleep(0.5)
end
=#

#=
import PyPlot: plt

x = [1,2,3]
y = [1,2,3]

plt.plot(x,y)
plt.show()
=#