module Wilsonloops
    export make_staple,Gaugeline, make_staple_and_loop
    using LaTeXStrings
    using LinearAlgebra
    import Base

    abstract type Gaugelink end

    struct GLink{Dim} <: Gaugelink
        direction::Int8
        position::NTuple{Dim,Int64}
    end

    struct Adjoint_GLink{Dim} <: Gaugelink
        parent::GLink{Dim}
    end

    function LinearAlgebra.adjoint(glink::GLink{Dim}) where Dim
        return Adjoint_GLink{Dim}(glink)
    end

    function LinearAlgebra.adjoint(glink::Adjoint_GLink{Dim}) where Dim
        return glink.parent
    end

    function get_direction(glink::GLink)
        return glink.direction
    end

    function get_direction(glink::Adjoint_GLink)
        return get_direction(glink.parent)
    end

    function get_position(glink::GLink)
        return glink.position
    end

    function get_position(glink::Adjoint_GLink)
        return get_position(glink.parent)
    end

    function set_position(glink::GLink{Dim},position) where Dim
        return GLink{Dim}(glink.direction,position)
    end

    function set_position(glink::Adjoint_GLink{Dim},position) where Dim
        return GLink{Dim}(glink.parent.direction,position)'
    end



    mutable struct Gaugeline{Dim}
        glinks::Array{Union{GLink{Dim},Adjoint_GLink{Dim}},1}
        Gaugeline(;Dim=4) = new{Dim}([])
        Gaugeline(glinks;Dim=4) = new{Dim}(glinks)
        function Gaugeline(segments_in::Array{Tuple{T,T},1};Dim=4) where T<: Integer
            segments = make_links(segments_in)
            numline = length(segments)
            glinks = Array{Union{GLink{Dim},Adjoint_GLink{Dim}},1}(undef,numline)
            position = zeros(Int64,Dim)
            for (i,segment) in enumerate(segments)
                dimension = segment[1]
                hoppingdirection = segment[2]

                if hoppingdirection == 1
                    glinks[i] =  GLink{Dim}(dimension ,Tuple(position))
                    position[dimension] += 1
                elseif hoppingdirection == -1
                    position[dimension] += -1
                    glinks[i] =  GLink{Dim}(dimension ,Tuple(position))'
                else
                    error("hoppingdirection in segment should be 1 or -1. But now $hoppingdirection")
                end
            end
            return new{Dim}(glinks)
        end
    end

    function Base.push!(w::Gaugeline,link)
        push!(w.glinks,link)
    end

    function Base.length(w::Gaugeline)
        return length(w.glinks)
    end

    function Base.getindex(w::Gaugeline,i)
        return w.glinks[i]
    end

    function Base.lastindex(w::Gaugeline)
        return length(w)
    end

    function LinearAlgebra.adjoint(w::Gaugeline{Dim}) where Dim
        wa = Gaugeline(;Dim=Dim)
        numlinks = length(w)
        for i=numlinks:-1:1
            push!(wa,w[i]')
        end
        return wa
    end

    function Base.display(ws::Array{<: Gaugeline{Dim},1}) where Dim
        for i=1:length(ws)
            println("$i-th loop")
            display(ws[i])
        end
    end

    function Base.display(w::Gaugeline{Dim}) where Dim
        outputstring = ""
        for (i,glink) in enumerate(w.glinks)
            direction  = get_direction(glink)
            nstring = "$(direction),n"

            dagornot = ifelse(typeof(glink) <: GLink,"","^{\\dagger}")
            position = get_position(glink)

            for μ=1:Dim
                m = position[μ]
                if m != 0
                    if abs(m)==1
                        if m >0
                            nstring = nstring*"+e_{$(μ)}"
                        else
                            nstring = nstring*"-e_{$(μ)}"
                        end
                    else
                        if m >0
                            nstring = nstring*"+$(m)e_{$(μ)}"
                        else
                            nstring = nstring*"-$(abs(m))e_{$(μ)}"
                        end
                    end
                end
            end
            outputstring = outputstring*"U_{$(nstring)}$(dagornot) "
        end
        println(outputstring)
        return outputstring
    end

    function make_staple(w::Gaugeline{Dim},μ) where Dim
        numlinks = length(w)
        linkindices = check_link(w,μ)
        numstaples = length(linkindices)
        staple = Array{typeof(w),1}(undef,numstaples)
        count = 0
        for (i,ith) in enumerate(linkindices)
            wi =Gaugeline(Dim=Dim)
            position = zero(collect(w[ith].position))
            position[w[ith].direction] += 1

            for j=ith+1:numlinks
                link = w[j]
                if typeof(link) <: GLink 
                    link_rev = set_position(link,Tuple(position))
                    position[get_direction(link)] += 1
                else
                    position[get_direction(link)] += -1
                    link_rev = set_position(link,Tuple(position))
                end
                push!(wi,link_rev)
            end

            for j=1:ith-1
                link = w[j]
                if typeof(link) <: GLink
                    link_rev = set_position(link,Tuple(position))
                    position[get_direction(link)] += 1
                else
                    position[get_dimension(link)] += -1
                    link_rev = set_direction(link,Tuple(position))
                end
                push!(wi,link_rev)
            end
            staple[i] = wi
            #println("μ = ",μ)
            #display(wi)
        end
        return staple

    end

    function make_staple_and_loop(w::Gaugeline{Dim},μ) where Dim
        C = make_staple(w,μ)
        append!(C,make_staple(w',μ))
        numstaple = length(C)
        CUdag = Array{typeof(w),1}(undef,numstaple)
        Udag = GLink{Dim}(μ,(0,0,0,0))'
        for i=1:numstaple
            CUdag[i] = deepcopy(C[i])'
            push!(CUdag[i],Udag)
            #CUdag[i] = CUdag[i]'
        end
        return CUdag
    end

    function check_link(w,μ)
        numlinks = length(w)
        linkindices=Int64[]
        for i=1:numlinks
            link = w[i]
            if typeof(link) <: GLink
                if link.direction == μ
                    append!(linkindices,i)
                end
            end
        end
        return linkindices
    end




    function make_links(segments::Array{Tuple{T,T},1}) where T <: Integer
        links = Tuple{Int8,Int8}[]
        for segment in segments
            s = sign(segment[2])
            if segment[2] == 0
                push!(links,(segment[1],0))
            else
            #@assert segment[2] != 0
                for i=1:abs(segment[2])
                    push!(links,(segment[1],s*1))
                end
            end
        end
        return links
    end

end