module Wilsonloops

    mutable struct Wilson_loop
        wilson_loop::Array{Tuple{Int8,Int8},1}
        origin::NTuple{4,Int8}
        Wilson_loop() = new([],(0,0,0,0))
        Wilson_loop(loop::Array{Tuple{Int,Int},1},origin) = new(loop,origin)
        Wilson_loop(loop,origin) = new([Tuple(loop)],origin)
        Wilson_loop(loop::Array{Tuple{Int,Int},1}) = new(loop,(0,0,0,0))
        Wilson_loop(loop) = new([Tuple(loop)],(0,0,0,0))
    end 

    function Base.display(w::Wilson_loop)
        println("origin = ",w.origin)
        print("[")
        for i=1:length(w)
            print(w.wilson_loop[i],",")
        end
        println("]\t")
    end



    function set_origin!(w::Wilson_loop,origin)
        w.origin = origin
    end


    mutable struct Wilson_loop_set
        #loopset::Set{Wilson_loop}
        loopset::Array{Wilson_loop,1}
        Wilson_loop_set() = new([])
    end

    function Base.display(w::Wilson_loop_set)
        for i=1:length(w)
            println("$i-th loops")
            display(w[i])
        end
        println("\t")
    end

    function Base.length(x::Wilson_loop_set)
        return length(x.loopset)
    end

    function Base.length(x::Wilson_loop)
        return length(x.wilson_loop)
    end

    function Base.getindex(g::Wilson_loop_set,i)
        return g.loopset[i]
    end

    function Base.getindex(g::Wilson_loop,i)
        return g.wilson_loop[i]
    end

    function Base.push!(g::Wilson_loop,link)
        push!(g.wilson_loop,Tuple(link))
        return
    end

    function Base.push!(g::Wilson_loop_set,w::Wilson_loop)
        push!(g.loopset,w)
        return
    end

    function make_loopforactions(couplinglist)
        loops = Array{Wilson_loop_set,1}(undef,length(couplinglist))
        for (i,name) in enumerate(couplinglist)
            if  name == "plaq"
                loops[i] = make_plaq()
            elseif name == "rect"
                loops[i] = make_rect()
            else
                error("$name is not supported!")
            end
        end
        return loops
    end


    function make_plaq()
        loops = Wilson_loop_set()
        origin = zeros(Int8,4)
        for μ=1:4
            for ν=1:4
                if ν == μ
                    continue
                end
                loop = make_links([(μ,1),(ν,1),(μ,-1),(ν,-1)])
                loop1 = Wilson_loop(loop,Tuple(origin))
                push!(loops,loop1)
            end
        end
        
        return loops
    end

    function make_plaq(origin)
        loops = Wilson_loop_set()
        for μ=1:4
            for ν=1:4
                if ν == μ
                    continue
                end
                loop = make_links([(μ,1),(ν,1),(μ,-1),(ν,-1)])
                loop1 = Wilson_loop(loop,Tuple(origin))
                push!(loops,loop1)
            end
        end
        
        return loops
    end

    function make_plaq(μ,ν)
        loops = Wilson_loop_set()
        origin = zeros(Int8,4)

        loop = make_links([(μ,1),(ν,1),(μ,-1),(ν,-1)])
        loop1 = Wilson_loop(loop,Tuple(origin))
        push!(loops,loop1)

        
        return loops
    end

    function make_plaq(μ,ν,origin)
        loops = Wilson_loop_set()
        #origin = zeros(Int8,4)

        loop = make_links([(μ,1),(ν,1),(μ,-1),(ν,-1)])
        loop1 = Wilson_loop(loop,Tuple(origin))
        push!(loops,loop1)

        
        return loops
    end

    function make_rect()




        loops = Wilson_loop_set()
        origin = zeros(Int8,4)

            

        for μ=1:4
            for ν=1:4
                if ν == μ
                    continue
                end
                #loop = make_links([(μ,1),(ν,1),(ν,1),(μ,-1),(ν,-1),(ν,-1)])
                loop = make_links([(μ,1),(ν,2),(μ,-1),(ν,-2)])
                #loop = make_links([(μ,2),(ν,1),(μ,-2),(ν,-1)])
                loop1 = Wilson_loop(loop,Tuple(origin))
                push!(loops,loop1)
            end
        end
        
        return loops
    end

    function make_plaq_staple(μ)
        loops = Wilson_loop_set()
        origin = zeros(Int8,4)
        origin[μ] = 1
        for ν=1:4
            if ν == μ
                continue
            end
            loop = make_links([(ν,1),(μ,-1),(ν,-1)])
            loop1 = Wilson_loop(loop,Tuple(origin))
            push!(loops,loop1)

            loop = make_links([(ν,-1),(μ,-1),(ν,1)])
            loop2 = Wilson_loop(loop,Tuple(origin))
            push!(loops,loop2)
        end
        return loops
    end

    function make_plaq_staple_prime(μ)
        loops = Wilson_loop_set()
        origin = zeros(Int8,4)
        #origin[μ] = 1
        for ν=1:4
            if ν == μ
                continue
            end
            loop = make_links([(ν,1),(μ,1),(ν,-1)])
            loop1 = Wilson_loop(loop,Tuple(origin))
            push!(loops,loop1)

            loop2 = Wilson_loop([(ν,-1),(μ,1),(ν,1)],Tuple(origin))
            push!(loops,loop2)
        end
        return loops
    end


    function calc_coordinate(w::Wilson_loop)
        numloops = length(w)
        coordinates = Array{NTuple{4,Int8},1}(undef,numloops)
        
        for k=1:numloops-1
            if k == 1
                coordinates[k] = w.origin
            end
            loopk = w[k]
            if loopk[1] == 1
                coordinates[k+1] = (coordinates[k][1]+loopk[2],coordinates[k][2],coordinates[k][3],coordinates[k][4])
            elseif loopk[1] == 2
                coordinates[k+1] = (coordinates[k][1],coordinates[k][2]+loopk[2],coordinates[k][3],coordinates[k][4])
            elseif loopk[1] == 3
                coordinates[k+1] = (coordinates[k][1],coordinates[k][2],coordinates[k][3]+loopk[2],coordinates[k][4])
            elseif loopk[1] == 4
                coordinates[k+1] = (coordinates[k][1],coordinates[k][2],coordinates[k][3],coordinates[k][4]+loopk[2])
            end
        end
        return coordinates
    end

    function calc_shift(w::Wilson_loop)
        numloops = length(w)
        coordinates = calc_coordinate(w)
        shifts = Array{NTuple{4,Int8},1}(undef,numloops)
        for k=1:numloops
            loopk = w[k]
            if loopk[2] == 1
                shifts[k] = coordinates[k]
            elseif loopk[2] == -1
                if loopk[1] == 1
                    shifts[k] = (coordinates[k][1]+loopk[2],coordinates[k][2],coordinates[k][3],coordinates[k][4])
                elseif loopk[1] == 2
                    shifts[k] = (coordinates[k][1],coordinates[k][2]+loopk[2],coordinates[k][3],coordinates[k][4])
                elseif loopk[1] == 3
                    shifts[k] = (coordinates[k][1],coordinates[k][2],coordinates[k][3]+loopk[2],coordinates[k][4])
                elseif loopk[1] == 4
                    shifts[k] = (coordinates[k][1],coordinates[k][2],coordinates[k][3],coordinates[k][4]+loopk[2])
                end
            end
        end
        return shifts
    end

    function make_links(segments::Array{Tuple{Int,Int},1})
        links = Tuple{Int,Int}[]
        for segment in segments
            s = sign(segment[2])
            @assert segment[2] != 0
            for i=1:abs(segment[2])
                push!(links,(segment[1],s*1))
            end
        end
        return links
    end

    function make_staples(w::Wilson_loop_set)
        num = length(w)
        staplesset = Array{Wilson_loop_set,1}(undef,4)
        tempset = Array{Array{Any,1},1}(undef,4)
        #tempset = Array{Set,1}(undef,4)
        for mu=1:4
            staplesset[mu] = Wilson_loop_set()
            #tempset[mu] = Set()
            tempset[mu] = []
        end

        
        for i=1:num
            staples = make_staples(w[i])
            for mu=1:4
                origin = zeros(Int8,4)
                origin[mu] = 1
                
                
                for k=1:length(staples[mu])
                    push!(tempset[mu],staples[mu][k])
                    #Wilson_loop(staples[mu][k],Tuple(origin)))
                end


                

                #for k=1:length(staples[mu])
                #    push!(staplesset[mu],Wilson_loop(staples[mu][k],Tuple(origin))  )
                #end
            end
        end

        for mu=1:4
            #display(tempset[mu])
            #exit()
            origin = zeros(Int8,4)
            origin[mu] = 1
            for temp in tempset[mu]
                push!(staplesset[mu],Wilson_loop(temp,Tuple(origin)))
            end
        end
        
        return staplesset
    end

    function make_staples(w::Wilson_loop)
        staples = make_staples(w.wilson_loop)
        return staples[1],staples[2],staples[3],staples[4]
    end

    function make_staples(links)
        numlinks = length(links)
        #println("num. of links: $numlinks")
        staples = Array{Array{Array{Tuple{Int,Int},1}}}(undef,4)
        for mu=1:4
            linkmu = []
            for (i,link) in enumerate(links)
                if link[1] == mu
                    append!(linkmu,i)
                end
            end
            #println(mu,"\t",linkmu)

            staple_mu = []

            if length(linkmu) == 0
                staple = Tuple{Int,Int}[]
                staple_mu = Array{Tuple{Int,Int},1}[]
            end

            staple_mu = Array{Array{Tuple{Int,Int},1}}(undef,length(linkmu))
            count = 0
            for i in linkmu
                count += 1
                staple = Array{Tuple{Int,Int}}(undef,numlinks-1)
                if i != numlinks
                    staple[1:numlinks-i] = links[i+1:end]
                end

                if i != 1
                    staple[numlinks-i+1:end] = links[1:i-1]
                end

                if links[i][2] == -1
                    staple = staple[end:-1:1]
                    for j = 1:numlinks-1
                        staple[j] = (staple[j][1],-staple[j][2])
                    end
                end             
                #println(staple)
                staple_mu[count] = staple[:]
            end
            staples[mu] = staple_mu
            #println(staple_mu)
        end

        #println(staples)
        return staples
    end




end