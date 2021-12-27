module Wilsonloops
    export make_staple,Gaugeline, make_staple_and_loop,derive_U,make_Cμ,
            make_plaq_staple,make_plaq, loops_staple_prime,get_position
    using LaTeXStrings
    using LinearAlgebra
    import Base

    abstract type Gaugelink{Dim} end

    struct GLink{Dim} <: Gaugelink{Dim}
        direction::Int8
        position::NTuple{Dim,Int64}
    end

    struct Adjoint_GLink{Dim} <: Gaugelink{Dim}
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

    struct DwDU{Dim,μ}
        parent::Gaugeline{Dim}
        insertindex::Int64
        position::NTuple{Dim,Int64}
        leftlinks::Gaugeline{Dim}
        rightlinks::Gaugeline{Dim}
    end

    function get_leftlinks(dw::DwDU)
        return dw.leftlinks
    end

    function get_rightlinks(dw::DwDU)
        return dw.rightlinks
    end



    function Base.push!(w::Gaugeline,link)
        push!(w.glinks,link)
    end

    function Base.append!(w::Gaugeline,a::Gaugeline)
        append!(w.glinks,a.glinks)
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

    function LinearAlgebra.adjoint(ws::Array{<: Gaugeline{Dim},1}) where Dim
        num = length(ws)
        wad = Array{eltype(ws),1}(undef,num)
        for i=1:num
            wad[i] = ws[i]'
        end
        return wad
    end


    function Base.show(io::IO,ws::Array{<: Gaugeline{Dim},1}) where Dim
        for i=1:length(ws)
            if i==1
                st ="st"
            elseif i==2
                st ="nd"
            elseif i==3
                st ="rd"
            else
                st ="th"
            end
            println(io,"$i-$st loop")
            show(io,ws[i])
            #display(io,ws[i])
        end
    end

    function Base.display(ws::Array{<: Gaugeline{Dim},1}) where Dim
        for i=1:length(ws)
            if i==1
                st ="st"
            elseif i==2
                st ="nd"
            elseif i==3
                st ="rd"
            else
                st ="th"
            end
            println("$i-$st loop")
            display(ws[i])
        end
    end

    function get_printstring_direction(glink::Gaugelink{Dim}) where Dim
        
        nstring = "n"
        
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
        return nstring
    end

    function get_printstring(glink::Gaugelink{Dim}) where Dim
        direction  = get_direction(glink)
        dagornot = ifelse(typeof(glink) <: GLink,"","^{\\dagger}")
        nstring = "$(direction),"*get_printstring_direction(glink)
        return "U_{$(nstring)}$(dagornot) "
    end

    function Base.show(io::IO,glink::Gaugelink{Dim}) where Dim
        outputstring =get_printstring(glink)
        println(io,outputstring)
    end

    function Base.display(w::Gaugeline{Dim}) where Dim
        outputstring = ""
        for (i,glink) in enumerate(w.glinks)
            outputstring = outputstring*get_printstring(glink)
        end
        println(outputstring)
        return outputstring
    end

    function Base.show(io::IO,w::Gaugeline{Dim}) where Dim
        outputstring = ""
        for (i,glink) in enumerate(w.glinks)
            outputstring = outputstring*get_printstring(glink)
        end
        println(io,outputstring)
        return outputstring
    end

    function make_staple(w::Gaugeline{Dim},μ) where Dim
        dwdUs = derive_U(w,μ)
        numstaples = length(dwdUs)
        staple = Array{typeof(w),1}(undef,numstaples)
        for i=1:numstaples
            wi =Gaugeline(Dim=Dim)
            append!(wi,get_rightlinks(dwdUs[i]))
            append!(wi,get_leftlinks(dwdUs[i]))
            staple[i] = wi
        end
        return staple
    end

    function make_Cμ(w::Gaugeline{Dim},μ) where Dim
        V1 = make_staple(w,μ)
        V2 = make_staple(w',μ)
        C = eltype(V1)[]
        for i=1:length(V1)
            push!(C,V1[i]')
        end
        for i=1:length(V2)
            push!(C,V2[i]')
        end
        return C
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

    """
        like U U U U -> U U otimes U 
    """
    function derive_U(w::Gaugeline{Dim},μ) where Dim
        numlinks = length(w)
        linkindices = check_link(w,μ)
        numstaples = length(linkindices)
        dwdU = Array{DwDU{Dim,μ},1}(undef,numstaples)

        for (i,ith) in enumerate(linkindices)
            #wi =Gaugeline(Dim=Dim)
            rightlinks = Gaugeline(Dim=Dim)
            leftlinks = Gaugeline(Dim=Dim)
            origin = w[ith].position
            position = zero(collect(origin))
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
                push!(rightlinks,link_rev)
            end

            for j=1:ith-1
                link = w[j]
                position = collect(get_position(link)) .- origin
                link_rev =  set_position(link,Tuple(position))
                push!(leftlinks,link_rev)
            end
            dwdU[i] = DwDU{Dim,μ}(w,ith,origin,leftlinks,rightlinks)
            #println("μ = ",μ)
            #display(wi)
        end
        return dwdU
    end

    function Base.display(dwdU::DwDU{Dim,μ}) where {Dim,μ}
        outputstring = ""
        if length(dwdU.leftlinks.glinks) == 0
            outputstring =outputstring*"I "
        else
            for glink in dwdU.leftlinks.glinks
                outputstring =outputstring*get_printstring(glink)
            end
        end

        outputstring = outputstring*" \\otimes "

        if length(dwdU.rightlinks.glinks) == 0
            outputstring =outputstring*"I "
        else
            for glink in dwdU.rightlinks.glinks
                outputstring =outputstring*get_printstring(glink)
            end
        end

        nstring = get_printstring_direction(dwdU.parent.glinks[dwdU.insertindex])

        outputstring = outputstring*"\\delta_{n,$(nstring)}"

        println(outputstring)
        return outputstring
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

    function make_plaq(μ,ν;Dim=4)
        return Gaugeline([(μ,1),(ν,1),(μ,-1),(ν,-1)],Dim = Dim)
    end

    function make_plaq(;Dim=4)
        loops = Gaugeline{Dim}[]
        for μ=1:Dim
            #for ν=1:4
            for ν=μ:Dim
                if ν == μ
                    continue
                end

                plaq = make_plaq(μ,ν,Dim=Dim)
                push!(loops,plaq)
            end
        end
        return loops
    end

    function make_rect(;Dim=4)
        loops = Gaugeline{Dim}[]
        for μ=1:4
            for ν=μ:4
                if ν == μ
                    continue
                end
                #loop = make_links([(μ,1),(ν,2),(μ,-1),(ν,-2)])

                loop1 = Gaugeline([(μ,1),(ν,2),(μ,-1),(ν,-2)],Dim = Dim)
                #loop1 = Wilson_loop([(μ,1),(ν,2),(μ,-1),(ν,-2)])
                #loop1 = Wilson_loop(loop,Tuple(origin))
                push!(loops,loop1)
                loop1 = Gaugeline([(μ,2),(ν,1),(μ,-2),(ν,-1)],Dim = Dim)
                #loop1 = Wilson_loop([(μ,2),(ν,1),(μ,-2),(ν,-1)])
                push!(loops,loop1)
            end
        end
        return loops
    end


    function make_plaq_staple(μ;Dim=4)
        loops = Gaugeline{Dim}[]
        plaqs = make_plaq(Dim = Dim)
        numplaq = length(plaqs)
        for i=1:numplaq
            plaq = plaqs[i]
            staples = make_staple(plaq,μ)
            for j=1:length(staples)
                push!(loops,staples[j])
            end

            plaqdag = plaqs[i]'
            staples = make_staple(plaqdag,μ)
            for j=1:length(staples)
                push!(loops,staples[j])
            end
        end

        return loops
    end

    function construct_staple_prime()
        loops_staple_prime = Dict{Tuple{Int8,Int8},Any}()
        for Dim=1:4
            for μ=1:Dim
                loops_staple_prime[(Dim,μ)] = make_plaq_staple(μ,Dim=Dim)'
            end
        end
        return loops_staple_prime
    end

    function construct_staple()
        loops_staple = Dict{Tuple{Int8,Int8},Any}()
        for Dim=1:4
            for μ=1:Dim
                loops_staple[(Dim,μ)] = make_plaq_staple(μ,Dim=Dim)
            end
        end
        return loops_staple
    end
    const loops_staple_prime = construct_staple_prime()
    const loops_staple = construct_staple()

end