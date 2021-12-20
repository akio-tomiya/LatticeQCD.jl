module FermionAction_module
    import ..Gaugefield:construct_smearing,Wilson_loop_set 
    import ..Gaugefield:Generator
    import ..Rhmc:RHMC

    abstract type FermiActionParam end

    Base.@kwdef struct FermiActionParam_Wilson{SmearingParam} <: FermiActionParam
        hop::Float64  = 0.141139#Hopping parameter
        r::Float64  = 1#Wilson term
        #eps::Float64 = 1e-8
        eps::Float64 = 1e-19
        Dirac_operator::String = "Wilson"
        MaxCGstep::Int64 = 3000 
        quench::Bool = false
        smearing::SmearingParam 

        function FermiActionParam_Wilson(hop,r,eps,Dirac_operator,MaxCGstep,quench
            ;smearingparameters = "nothing",
            loops_list = nothing,
            coefficients  = nothing,
            numlayers = 1,
            L = nothing)

            smearing = construct_smearing(smearingparameters,loops_list,L,coefficients,numlayers)
            #smearing = Nosmearing()
            #if smearingparameters == nothing
            #    smearing = Nosmearing()
            #end
            return new{typeof(smearing)}(hop,r,eps,Dirac_operator,MaxCGstep,quench,smearing)
        end

    end

    const Clover_coefficient  = 1.5612

    Base.@kwdef struct FermiActionParam_WilsonClover{SmearingParam} <: FermiActionParam
        hop::Float64  = 0.141139#Hopping parameter
        r::Float64  = 1#Wilson term
        #eps::Float64 = 1e-8
        eps::Float64 = 1e-19
        Dirac_operator::String = "WilsonClover"
        MaxCGstep::Int64 = 3000 
        Clover_coefficient::Float64 = Clover_coefficient 
        #CloverFμν::Array{ComplexF64,4}
        internal_flags::Array{Bool,1}
        inn_table::Array{Int64,3}
        _ftmp_vectors::Array{Array{ComplexF64,3},1}
        _is1::Array{Int64,1}
        _is2::Array{Int64,1}
        quench::Bool = false
        SUNgenerator::Union{Nothing,Generator}
        _cloverloops::Array{Wilson_loop_set,2}
        smearing::SmearingParam 

        function FermiActionParam_WilsonClover(hop,r,eps,Dirac_operator,MaxCGstep,
            Clover_coefficient,internal_flags,inn_table,_ftmp_vectors,
            _is1,_is2,
            quench, SUNgenerator,_cloverloops;smearingparameters = "nothing",
            loops_list = nothing,
            coefficients  = nothing,
            numlayers = 1,
            L = nothing)

            smearing = construct_smearing(smearingparameters,loops_list,L,coefficients,numlayers)
            
            #if smearingparameters == nothing
            #    smearing = Nosmearing()
            #end
            #smearing = Nosmearing()

            return new{typeof(smearing)}(hop,r,eps,Dirac_operator,MaxCGstep,
            Clover_coefficient,internal_flags,inn_table,_ftmp_vectors,
            _is1,_is2,
            quench, SUNgenerator,_cloverloops,smearing)
        end
    end


    

    Base.@kwdef struct FermiActionParam_Staggered{SmearingParam} <: FermiActionParam
        mass::Float64 = 0.5
        eps::Float64 = 1e-19
        Dirac_operator::String = "Staggered"
        MaxCGstep::Int64 = 3000 
        quench::Bool = false
        Nf::Int8 = 4
        rhmc_action::Union{Nothing,RHMC}
        rhmc_MD::Union{Nothing,RHMC}
        smearing::SmearingParam



        function FermiActionParam_Staggered(
            mass,
            eps,
            Dirac_operator,
            MaxCGstep,
            quench,
            Nf;smearingparameters = "nothing",
            loops_list = nothing,
            coefficients  = nothing,
            numlayers = 1,
            L = nothing
            ) where T <: Real

            if Nf == 4 || Nf == 8 # 8 flavors if phi (Mdag M)^{-1} phi
                rhmc_action = nothing
                rhmc_MD = nothing
            else
                #for action: r_action
                #Nf = 8 -> alpha = 1 -> power x^{1/2} 
                #Nf = 2 -> alpha = 1/4 -> power x^1/8 
                #Nf = 1 -> alpha = 1/8  -> power x^1/16 
                order = Nf //16

                rhmc_action = RHMC(order,n=15)

                #for MD: r_MD
                #Nf = 8 -> alpha = 1 -> power x^{1} 
                #Nf = 2 -> alpha = 1/4 -> power x^1/4 
                #Nf = 1 -> alpha = 1/8  -> power x^1/8 
                order = Nf // 8
                #rhmcorder = 8 ÷ Nf
                rhmc_MD = RHMC(order,n=10)
            end

            smearing = construct_smearing(smearingparameters,loops_list,L,coefficients,numlayers)


            return new{typeof(smearing)}(
            mass,
            eps,
            Dirac_operator,
            MaxCGstep,
            quench,
            Nf,
            rhmc_action,
            rhmc_MD,
            smearing
            )
        end
    end

end