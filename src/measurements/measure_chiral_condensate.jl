module Measure_chiral_condensate_module
    using LinearAlgebra
    import ..AbstractMeasurement_module:AbstractMeasurement,measure,set_params,
            set_Wilson,set_WilsonClover,set_Staggered
    import ..Gaugefield:AbstractGaugefields,
                                        Traceless_antihermitian!,Traceless_antihermitian,
                                        initialize_TA_Gaugefields
    import ..Gaugefields:calculate_Plaquette,calculate_Polyakov_loop,substitute_U!,calculate_energy_density
    import ..Verbose_print:Verbose_level,Verbose_3,Verbose_2,Verbose_1,println_verbose3,println_verbose2,println_verbose1,
            print_verbose1,print_verbose2,print_verbose3
    import ..Gaugefield:Wilson_loop,Wilson_loop_set,calc_loopset_μν_name            
    import ..Gaugefield:Loops,evaluate_loops,TA_Gaugefields
    import ..Smearing:gradientflow!,calc_stout!,calc_fatlink_APE!,calc_stout,calc_fatlink_APE,calc_multihit!
    import ..Actions:GaugeActionParam_autogenerator,GaugeActionParam
    import ..Actions:FermiActionParam_Wilson,FermiActionParam_Staggered,FermiActionParam_WilsonClover,
                FermiActionParam




    mutable struct Measure_chiral_condensate{T} <: AbstractMeasurement
        filename::String
        fp::IOStream
        tempU::Array{T,1}
        printvalues::Bool
        fparam::FermiActionParam

        function Measure_chiral_condensate(filename,
                    U::Array{T,1},params;printvalues = true) where T
            fp = open(filename,"w")
            
            tempU = Array{T,1}(undef,3)
            for i=1:3
                tempU[i] = similar(U[1])
            end

            if haskey(params,"fermiontype")
                fermiontype = params["fermiontype"] 
            else
                fermiontype = nothing
            end

            if haskey(params,"eps")
                eps = params["eps"]
            else
                eps = 1e-16
            end

            if haskey(params,"MaxCGstep")
                MaxCGstep = params["MaxCGstep"]
            else
                MaxCGstep = 5000
            end

            if fermiontype  == "Wilson"
                fparam = set_Wilson(params)
            elseif fermiontype == "WilsonClover"
                fparam = set_WilsonClover(params)
            elseif fermiontype == "Staggered"
                fparam = set_Staggered(params)
            elseif fermiontype == nothing
                fparam = nothing
            else
                error("$fermiontype is not supported. use Wilson or Staggered")
            end


            m = new{T}(filename,fp,tempU,printvalues,fparam
            )
            finalizer(m) do m
                close(m.fp)
            end
            return m
        end

    end


    function measure(m::M,itrj,U::Array{<: AbstractGaugefields{NC,Dim},1};verbose = Verbose_2()) where {M <: Measure_chiral_condensate,NC,Dim}

        error("not implemented")
        return 
    end

end