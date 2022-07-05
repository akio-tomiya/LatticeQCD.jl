const BoundaryCondition_4D_default = [1, 1, 1, -1]
const BoundaryCondition_2D_default = [1, -1]

mutable struct Chiral_condensate_measurement{Dim,TG,TD,TF,TF_vec} <: AbstractMeasurement
    filename::String
    _temporary_gaugefields::Vector{TG}
    Dim::Int8
    #factor::Float64
    verbose_print::Union{Verbose_print,Nothing}
    printvalues::Bool
    D::TD
    fermi_action::TF
    _temporary_fermionfields::Vector{TF_vec}
    Nr::Int64
    factor::Float64

    function Chiral_condensate_measurement(
        U::Vector{T};
        filename = nothing,
        verbose_level = 2,
        printvalues = true,
        fermiontype = "Staggered",
        mass = 0.1,
        Nf = 2,
        κ = 1,
        r = 1,
        L5 = 2,
        M = -1,
        eps_CG = 1e-14,
        MaxCGstep = 3000,
        BoundaryCondition = nothing,
        Nr = 10,
    ) where {T}

        Dim = length(U)
        if BoundaryCondition == nothing
            if Dim == 4
                boundarycondition = BoundaryCondition_4D_default
            elseif Dim == 2
                boundarycondition = BoundaryCondition_2D_default
            end
        else
            boundarycondition = BoundaryCondition
        end
        Nfbase = 1
        factor = 1
        params = Dict()
        parameters_action = Dict()
        if fermiontype == "Staggered"
            x = Initialize_pseudofermion_fields(U[1], "staggered")
            params["Dirac_operator"] = "staggered"
            params["mass"] = mass
            parameters_action["Nf"] = Nf
            Nfbase = 4
            #Nfbase = ifelse( m.fparam.Dirac_operator == "Staggered",4,1)
            factor = Nf / Nfbase

        elseif fermiontype == "Wilson"
            x = Initialize_pseudofermion_fields(U[1], "Wilson", nowing = true)
            params["Dirac_operator"] = "Wilson"
            params["κ"] = κ
            params["r"] = r
            params["faster version"] = true
        elseif fermiontype == "Domainwall"
            params["Dirac_operator"] = "Domainwall"
            params["mass"] = mass
            params["L5"] = L5
            params["M"] = M
            x = Initialize_pseudofermion_fields(U[1], "Domainwall", L5 = L5)
        else
            error(
                "fermion type $fermiontype is not supported in chiral condensate measurement",
            )
        end

        params["eps_CG"] = eps_CG
        params["verbose_level"] = verbose_level
        params["MaxCGstep"] = MaxCGstep
        params["boundarycondition"] = boundarycondition

        D = Dirac_operator(U, x, params)
        fermi_action = FermiAction(D, parameters_action)
        TD = typeof(D)
        TF = typeof(fermi_action)



        myrank = get_myrank(U)
        #=
        if U[1].mpi == false
            myrank = 0
        else
            myrank = U[1].myrank
        end
        =#
        if printvalues
            verbose_print = Verbose_print(verbose_level, myid = myrank, filename = filename)
        else
            verbose_print = nothing
        end



        numg = 1
        _temporary_gaugefields = Vector{T}(undef, numg)
        _temporary_gaugefields[1] = similar(U[1])
        for i = 2:numg
            _temporary_gaugefields[i] = similar(U[1])
        end

        numf = 2
        TF_vec = typeof(x)
        _temporary_fermionfields = Vector{TF_vec}(undef, numf)
        for i = 1:numf
            _temporary_fermionfields[i] = similar(x)
        end

        return new{Dim,T,TD,TF,TF_vec}(
            filename,
            _temporary_gaugefields,
            Dim,
            verbose_print,
            printvalues,
            D,#::TD
            fermi_action,#::TF,
            _temporary_fermionfields,
            Nr,
            factor,
        )

    end



end

function measure(
    m::M,
    itrj,
    U::Array{<:AbstractGaugefields{NC,Dim},1};
    additional_string = "",
) where {M<:Chiral_condensate_measurement,NC,Dim}
    temps_fermi = get_temporary_fermionfields(m)
    p = temps_fermi[1]
    r = temps_fermi[2]
    D = m.D(U)
    pbp = 0.0
    #Nr = 100
    Nr = m.Nr

    for ir = 1:Nr
        clear_fermion!(p)
        Z4_distribution_fermi!(r)
        solve_DinvX!(p, D, r)
        tmp = dot(r, p) # hermitian inner product

        if m.printvalues
            #println_verbose_level2(U[1],"# $itrj $ir $(real(tmp)/U[1].NV) # itrj irand chiralcond")
            println_verbose_level2(
                m.verbose_print,
                "# $itrj $ir $additional_string $(real(tmp)/U[1].NV) # itrj irand chiralcond",
            )
        end
        pbp += tmp
    end

    pbp_value = real(pbp / Nr) / U[1].NV * m.factor

    if m.printvalues
        println_verbose_level1(m.verbose_print, "$itrj $pbp_value # pbp Nr=$Nr")
        flush(stdout)
    end

    return pbp_value
end



#=
"""
c-------------------------------------------------c
c     Random number function Z4  Noise
c     https://arxiv.org/pdf/1611.01193.pdf
c-------------------------------------------------c
    """
    function Z4_distribution_fermion!(x::AbstractFermionfields_4D{NC})  where NC
        NX = x.NX
        NY = x.NY
        NZ = x.NZ
        NT = x.NT
        n6 = size(x.f)[6]
        θ = 0.0
        N::Int32 = 4
        Ninv = Float64(1/N)
        for ialpha = 1:n6
            for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        for ix=1:NX
                            @inbounds @simd for ic=1:NC
                                θ = Float64(rand(0:N-1))*π*Ninv # r \in [0,π/4,2π/4,3π/4]
                                x[ic,ix,iy,iz,it,ialpha] = cos(θ)+im*sin(θ) 
                            end
                        end
                    end
                end
            end
        end

        set_wing_fermion!(x)

        return
    end

    =#
