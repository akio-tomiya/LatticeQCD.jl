module Staple_generator

    mutable struct Wilson_loop_term
        wilson_loops::Array{Array{Tuple{Int,Int},1}}
        Wilson_loop_term() = new([])
    end 

    function Base.push!(w::Wilson_loop_term,loop)
        push!(w.wilson_loops,loop)
        return
    end

    mutable struct Generator
        wilson_loop_terms::Array{Wilson_loop_term,1}
        Generator() = new([])
    end


    function Base.getindex(g::Generator,i)
        return g.wilson_loop_terms[i]
    end

    function Base.getindex(w::Wilson_loop_term,i)
        return w.wilson_loops[i]
    end

    function Base.push!(g::Generator,w::Wilson_loop_term)
        push!(g.wilson_loop_terms,w)
        return
    end

    function Base.length(g::Generator)
        return length(g.wilson_loop_terms)
    end

    function Base.length(w::Wilson_loop_term)
        return length(w.wilson_loops)
    end

    function rects(len)
        rect = Wilson_loop_term() 
        #segments_group = Array{Tuple{Int,Int},1}[]
        for mu=1:4
            for nu=1:4
                if nu == mu
                    continue
                end
                a = [(mu,1),(nu,len),(mu,-1),(nu,-len)]
                push!(rect,a)
                #push!(segments_group,a)
            end
        end
        #return segments_group
        return rect
    end

    function crown()
        segments_group = Array{Tuple{Int,Int},1}[]
        nu1 = [1,2,3]
        nu2 = [1,2,4]
        nu3 = [2,3,4]
        nu4 = [1,3,4]

        numurho = [
        (nu1[1],nu1[2],nu1[3]),(nu1[1],nu1[3],nu1[2]),(nu1[2],nu1[1],nu1[3]),(nu1[2],nu1[3],nu1[1]),(nu1[3],nu1[1],nu1[2]),(nu1[3],nu1[2],nu1[1]),
        (nu2[1],nu2[2],nu2[3]),(nu2[1],nu2[3],nu2[2]),(nu2[2],nu2[1],nu2[3]),(nu2[2],nu2[3],nu2[1]),(nu2[3],nu2[1],nu2[2]),(nu2[3],nu2[2],nu2[1]),
        (nu3[1],nu3[2],nu3[3]),(nu3[1],nu3[3],nu3[2]),(nu3[2],nu3[1],nu3[3]),(nu3[2],nu3[3],nu3[1]),(nu3[3],nu3[1],nu3[2]),(nu3[3],nu3[2],nu3[1]),
        (nu4[1],nu4[2],nu4[3]),(nu4[1],nu4[3],nu4[2]),(nu4[2],nu4[1],nu4[3]),(nu4[2],nu4[3],nu4[1]),(nu4[3],nu4[1],nu4[2]),(nu4[3],nu4[2],nu4[1])]
        segments_group = Array{Tuple{Int,Int},1}[]
        for indices in numurho
            nu,mu,rho = indices
            #println(nu,"\t",mu,"\t",rho)
            a = [(mu,1),(nu,1),(rho,1),(mu,-1),(nu,-1),(rho,-1)]
            push!(segments_group,a)
        end
        return segments_group
    end

    function chair()
        segments_group = Array{Tuple{Int,Int},1}[]
        nu1 = [1,2,3]
        nu2 = [1,2,4]
        nu3 = [2,3,4]
        nu4 = [1,3,4]

        numurho = [
        (nu1[1],nu1[2],nu1[3]),(nu1[1],nu1[3],nu1[2]),(nu1[2],nu1[1],nu1[3]),(nu1[2],nu1[3],nu1[1]),(nu1[3],nu1[1],nu1[2]),(nu1[3],nu1[2],nu1[1]),
        (nu2[1],nu2[2],nu2[3]),(nu2[1],nu2[3],nu2[2]),(nu2[2],nu2[1],nu2[3]),(nu2[2],nu2[3],nu2[1]),(nu2[3],nu2[1],nu2[2]),(nu2[3],nu2[2],nu2[1]),
        (nu3[1],nu3[2],nu3[3]),(nu3[1],nu3[3],nu3[2]),(nu3[2],nu3[1],nu3[3]),(nu3[2],nu3[3],nu3[1]),(nu3[3],nu3[1],nu3[2]),(nu3[3],nu3[2],nu3[1]),
        (nu4[1],nu4[2],nu4[3]),(nu4[1],nu4[3],nu4[2]),(nu4[2],nu4[1],nu4[3]),(nu4[2],nu4[3],nu4[1]),(nu4[3],nu4[1],nu4[2]),(nu4[3],nu4[2],nu4[1])]
        segments_group = Array{Tuple{Int,Int},1}[]
        for indices in numurho
            nu,mu,rho = indices
            #println(nu,"\t",mu,"\t",rho)
            a = [(mu,1),(nu,1),(rho,1),(mu,-1),(rho,-1),(nu,-1)]
            push!(segments_group,a)
        end
        return segments_group
    end


    function write_actions(fp,header,g::Generator)
        numactions = length(g)
        for i=1:numactions
            header_i = "$(i)th_"*header
            make_functions(fp,header_i,g[1])
        end
    end

    function make_functions(fp,header,w::Wilson_loop_term)
        numloops = length(w)
        for i=1:numloops
            segments = w[i]
            println(fp,"# Loop $i is ",segments)
            check_loop(segments)
            links = make_links(segments)
            println(fp,"# Loop is ", links)
            staples = make_staples(links)
            numstaples = make_function(fp,i,header,staples)
            #numstaples = make_function(fp,i,header,staples)
        end
    end

    function make_function(fp,i,header,staples)
        tabs = ""
        println(fp,"function calc_staple_"*header*"_num$(i)(U,mu,temps)")
        tabs = tabs*"\t"
        println(fp,tabs*"temp1 = temps[1]")
        println(fp,tabs*"temp2 = temps[2]")
        println(fp,tabs*"temp3 = temps[3]")
        println(fp,tabs*"temp4 = temps[4]")
        println(fp,tabs*"staple = temps[5]")

        nums = print_staples(fp,tabs,staples)
        println(fp,tabs*"return V")
        println(fp,"end")
        return nums
    end

    


    function make_ustring(link)
        if link[2] == 1
            return "U[nn,$(link[1])]"
        elseif link[2] == -1
            return "U[nn,$(link[1])]'"
        end
    end

    function check_loop(segments::Array{Tuple{Int,Int},1})
        position = zeros(Int64,4)
        for segment in segments
            position[segment[1]] += segment[2]
        end
        if sum(position) != 0
            println("Warning! This is not a closed loop! End point is at $position")
            println("Warning! $segments")
        end
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


    function default_wilson_loops()
        generator = Generator()
        len = 3
        rect = Staple_generator.rects(len)
        push!(generator,rect)
        return generator

        wilson_loops = Array{Array{Tuple{Int,Int},1}}[]
        len = 3
        wilson_loops = Array{Array{Tuple{Int,Int},1}}[]
        wilson_loopset1 = rects(len)
        wilson_loopset2 = crown()
        wilson_loopset3 = chair()
        push!(wilson_loops,wilson_loopset1)
        push!(wilson_loops,wilson_loopset2)
        push!(wilson_loops,wilson_loopset3)

        return segments_groups
    end




end

#=
using .Staple_generator
len = 3
println(typeof(Staple_generator.rects(len)))
println(Staple_generator.rects(len))

g = Staple_generator.default_wilson_loops()
fp = open("teststaple.jl","w")
header = "original"
Staple_generator.write_actions(fp,header,g)
close(fp)

=#