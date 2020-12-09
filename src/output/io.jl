module IOmodule
    using JLD
    import ..Gaugefields:GaugeFields

    function saveU(filename,x::Array{T,1}) where T <: GaugeFields
        NX=x[1].NX
        NY=x[1].NY
        NZ=x[1].NZ
        NT=x[1].NT
        NC=x[1].NC
        NDW=x[1].NDW
        NV=x[1].NV
        save(filename,"NX",NX,"NY",NY,"NZ",NZ,"NT",NT,"NC",NC,"NDW",NDW,"NV",NV,"U",x)
        return
    end

    function loadU!(filename,U) where T <: GaugeFields
        NX=load(filename, "NX")
        NY=load(filename, "NY")
        NZ=load(filename, "NZ")
        NT=load(filename, "NT")
        NC=load(filename, "NC")
        NDW=load(filename, "NDW")
        NV=load(filename, "NV")

        @assert NX==U[1].NX
        @assert NY==U[1].NY
        @assert NZ==U[1].NZ
        @assert NT==U[1].NT
        @assert NC==U[1].NC
        @assert NDW==U[1].NDW
        @assert NV==U[1].NV


        Unew = load(filename, "U")
        for μ=1:4
            U[μ] = Unew[μ]
        end
        return 
    end

    function loadU(filename) where T <: GaugeFields
        NX=load(filename, "NX")
        NY=load(filename, "NY")
        NZ=load(filename, "NZ")
        NT=load(filename, "NT")
        NC=load(filename, "NC")
        NDW=load(filename, "NDW")
        NV=load(filename, "NV")

        return  load(filename, "U")

    end

    function loadU(filename,NX,NY,NZ,NT,NC) where T <: GaugeFields
        @assert NX == load(filename, "NX") "NX in file $filename is $(load(filename, "NX")) but NX = $NX is set" 
        @assert NY == load(filename, "NY") "NY in file $filename is $(load(filename, "NY")) but NY = $NY is set" 
        @assert NZ == load(filename, "NZ") "NZ in file $filename is $(load(filename, "NZ")) but NZ = $NZ is set" 
        @assert NT == load(filename, "NT") "NT in file $filename is $(load(filename, "NY")) but NT = $NT is set" 
        @assert NC == load(filename, "NC") "NC in file $filename is $(load(filename, "NC")) but NC = $NC is set" 
        NDW=load(filename, "NDW")
        NV=load(filename, "NV")

        return  load(filename, "U")
    end


    abstract type IOFormat end
end