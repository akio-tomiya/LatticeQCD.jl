module Analyze
    #using Plots
    using Requires

    function analyze(dirname)
        return Analyze_data(dirname)
    end



    mutable struct Analyze_data
        plaquette::Array{Float64,1}
        polyakov::Array{ComplexF64,1}
        trjs::Array{Int64,1}
        numdata::Int64
        dirname::String

        function Analyze_data(dirname)
            plaquette = Float64[]
            polyakov = ComplexF64[]
            trjs = Int64[]
            filename_plaquette = pwd()*"/"*dirname*"/Plaquette.txt"
            filename_polyakov = pwd()*"/"*dirname*"/Polyakov_loop.txt"
            open(filename_plaquette) do f
                for i in eachline(f)
                    if occursin("plaq",i)
                        tmp = split(i)
                        itrj = parse(Int,tmp[1])
                        plaq = parse(Float64,tmp[2])
                        append!(trjs,itrj)
                        append!(plaquette,plaq)
                    end
                    
                end
            end

            open(filename_polyakov) do f
                for i in eachline(f)
                    if occursin("poly",i)
                        tmp = split(i)
                        itrj = parse(Int,tmp[1])
                        poly = parse(Float64,tmp[2])+im*parse(Float64,tmp[3])
                        append!(polyakov,poly)
                    end
                    
                end
            end

            return new(plaquette,polyakov,trjs,length(trjs),dirname)
        end
    end

    function get_plaquette(dirname::String)
        a = analyze(dirname)
        return a.plaquette
    end

    function get_polyakov(dirname::String)
        a = analyze(dirname)
        return a.polyakov
    end

    function get_trjs(dirname::String)
        a = analyze(dirname)
        return a.trjs
    end

    function get_plaquette(a::Analyze_data)
        return a.plaquette
    end

    function get_polyakov(a::Analyze_data)
        return a.polyakov
    end

    function get_trjs(a::Analyze_data)
        return a.trjs
    end

    function get_dirname(a::Analyze_data)
        return a.dirname
    end

    function __init__()
        @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin     
            using Plots        
            function plot_plaquette(dirname::String)
                a = analyze(dirname)
                plot_plaquette(a)
            end

            function plot_plaquette(a::Analyze_data)
                plt1 = histogram([0],label=nothing) #plot1
                plt2 = plot([], [],label=nothing) #plot2
                ylabel!("Plaquette")
                xlabel!("MC time")
                plot(plt1, plt2,  layout = 2)
                hist_plaq=[]

                Nsteps = a.numdata

                function plot_refresh!(plt1,plt2,hist_plaq,plaq,itrj)
                    bins=round(Int,log(itrj)*5+1)
                    append!(hist_plaq, plaq)

                    plt1 = histogram(hist_plaq,bins=bins,label=nothing) #plot1
                    xlabel!("Plaquette")

                    push!(plt2, 1, itrj, plaq)

                    plot(plt1, plt2, layout = 2)

                end

                for itrj=1:Nsteps
                    plaq = a.plaquette[itrj]

                    plot_refresh!(plt1,plt2,hist_plaq,plaq,itrj)

                    flush(stdout)
                end

                savefig(pwd()*"/"*a.dirname*"/Plaquette.pdf")
            end

            function plot_polyakov(dirname::String)
                a = analyze(dirname)
                plot_polyakov(a)
            end

            function plot_polyakov(a::Analyze_data)
                plt1 = histogram([0],label=nothing) #plot1
                plt2 = plot([], [],label=nothing) #plot2
                ylabel!("Polyakov loop")
                xlabel!("MC time")
                plot(plt1, plt2,  layout = 2)
                hist_poly=[]

                Nsteps = a.numdata

                function plot_refresh!(plt1,plt2,hist_poly,poly,itrj)
                    bins=round(Int,log(itrj)*5+1)
                    append!(hist_poly,poly)

                    plt1 = histogram(hist_poly,bins=bins,label=nothing) #plot1
                    xlabel!("Polyakov loop")

                    push!(plt2, 1, itrj, poly)

                    plot(plt1, plt2, layout = 2)

                end

                for itrj=1:Nsteps
                    poly = a.polyakov[itrj]
                    plot_refresh!(plt1,plt2,hist_poly,abs(poly),itrj)

                    flush(stdout)
                end
                savefig(pwd()*"/"*a.dirname*"/Polyakov.pdf")
            end

            function plot_plaq_and_poly(dirname::String)
                a = analyze(dirname)
                plot_plaq_and_poly(a)
            end

            function plot_plaq_and_poly(a::Analyze_data)
                plt1 = histogram([0],label=nothing) #plot1
                plt2 = plot([], [],label=nothing) #plot2
                ylabel!("Plaquette")
                xlabel!("MC time")
                plt3 = histogram([0],label=nothing) #plot3
                plt4 = plot([], [],label=nothing) #plot4
                ylabel!("Polyakov loop")
                xlabel!("MC time")
                plot(plt1, plt2, plt3, plt4, layout = 4)

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
                    #gui()
                end

                Nsteps = a.numdata

                for itrj=1:Nsteps
                    poly = a.polyakov[itrj]
                    plaq = a.plaquette[itrj]
                    plot_refresh!(plt1,plt2,plt3,plt4,hist_plaq,hist_poly, plaq, abs(poly),itrj)

                    flush(stdout)
                end
                savefig(pwd()*"/"*a.dirname*"/Plaquette_Polyakov.pdf")
            end
        end
    end

end