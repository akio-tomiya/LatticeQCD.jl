mutable struct Stoutsmearing <: Abstractsmearing
    staples_for_stout::Array{Array{Wilson_loop_set,1},1}
    tensor_derivative::Array{Tensor_derivative_set,1}
    staples_for_stout_dag::Array{Array{Wilson_loop_set,1},1}
    tensor_derivative_dag::Array{Tensor_derivative_set,1}
    ρs::Array{Array{Float64,1},1}
end

struct STOUT_dataset{Dim}
    closedloop::Vector{Wilsonline{Dim}}
    Cμ::Vector{Vector{Wilsonline{Dim}}}
    dCμdUν::Matrix{Vector{DwDU{Dim}}} #(mu,nu) num loops
    dCμdUνdag::Matrix{Vector{DwDU{Dim}}} #(mu,nu) num loops
    #dCmudUnu::Dict{Tuple{Int64,Int64,Int64},Array{DwDU{Dim},1}}
end

struct STOUT_Layer{Dim} <: CovLayer{Dim}
    ρs::Vector{Float64}
    dataset::Vector{STOUT_dataset{Dim}}
end

function Base.length(layer::STOUT_Layer)
    return length(layer.ρs)
end

function get_Cμ(layer::STOUT_Layer,i)
    return layer.dataset[i].Cμ
end

function get_dCμdUν(layer::STOUT_Layer,i)
    return layer.dataset[i].dCμdUν
end

function get_dCμdUνdag(layer::STOUT_Layer,i)
    return layer.dataset[i].dCμdUνdag
end

function get_ρ(layer::STOUT_Layer,i)
    return layer.ρs[i]
end




function CovNeuralnet_STOUT(loops_smearing,ρs,L;Dim=4)

    numlayers = length(ρs)
    layers = Array{CovLayer,1}(undef,numlayers)
    for i=1:numlayers
        stout_layer = STOUT_Layer(loops_smearing,ρs[i],L,Dim=Dim)
        layers[i] = stout_layer
    end

    return CovNeuralnet{Dim}(numlayers,layers)#,_temp_gaugefields,_temp_TA_gaugefields)

end

#=
closedloops -> sum_{nm} P_nm
Cmu = dP_nm/dUmu
dCmu/dUnu = (d/dUnu) dP_nm/dUmu
=#

function STOUT_dataset(closedloops;Dim=4)
    #Ci = #Dict{Tuple{Int64,Int64},Array{Wilsonline{Dim},1}}[]
    #dCmudUnu = #Dict{Tuple{Int64,Int64,Int64},Array{DwDU{Dim},1}}[]
    num = length(closedloops) #number of loops 

    Cμs = Vector{Vector{Wilsonline{Dim}}}(undef,Dim)
    for μ=1:Dim
        Cμs[μ] = Vector{Wilsonline{Dim}}[] #set of staples. In the case of plaq, there are six staples. 
    end

    for i=1:num
        glinks = closedloops[i]
        for μ=1:Dim
            Cμ = make_Cμ(glinks,μ)
            for j=1:length(Cμ)
                push!(Cμs[μ],Cμ[j])
            end
        end
    end

    CmudUnu = Matrix{Vector{DwDU{Dim}}}(undef,Dim,Dim)
    CmudUnudag = Matrix{Vector{DwDU{Dim}}}(undef,Dim,Dim)


    for ν=1:Dim
        for μ=1:Dim
            CmudUnu[μ,ν] = Vector{DwDU{Dim}}[] 
            CmudUnudag[μ,ν] = Vector{DwDU{Dim}}[] 
            Cμ = Cμs[μ]
            numCμ = length(Cμ)
            for j=1:numCμ
                Cμj = Cμ[j]
                dCμjν = derive_U(Cμj,ν)
                numdCμjν = length(dCμjν)
                for k=1:numdCμjν
                    push!(CmudUnu[μ,ν],dCμjν[k])
                end

                dCμjνdag = derive_Udag(Cμj,ν)
                numdCμjνdag = length(dCμjνdag)
                for k=1:numdCμjνdag
                    push!(CmudUnudag[μ,ν],dCμjνdag[k])
                end

            end
            #println("dC$(μ)/dU$(ν): ")
            #show(CmudUnu[μ,ν])
        end
    end


    for μ=1:Dim
        println("μ = $μ")
        show(Cμs[μ])
        for ν=1:Dim
            println("dC$(μ)/dU$(ν): ")
            show(CmudUnu[μ,ν])
        end

        for ν=1:Dim
            println("dC$(μ)/dUdag$(ν): ")
            show(CmudUnudag[μ,ν])
        end
    end


    return STOUT_dataset{Dim}(closedloops,Cμs,CmudUnu,CmudUnudag)
end

function STOUT_Layer(loops_smearing,ρ,L;Dim=4)
    num = length(loops_smearing)
    loopset = make_loopforactions(loops_smearing,L)
    dataset = Array{STOUT_dataset{Dim},1}(undef,num)
    for i=1:num
        closedloops = loopset[i] #one of loopset, like plaq. There are several loops. 
        dataset[i] =  STOUT_dataset(closedloops,Dim=Dim)
    end

    return STOUT_Layer{Dim}(ρ,dataset)
end

"""
δ_prev = δ_current*exp(Q) - C^+ Λ 
        + sum_{μ',m} sum_{i} [B^i_{μ,m} U_{μ'}^+ Λ_{μ',m} A_{μ,m} - bar{B}_{μ,m} Λ_{μ',m}U_{μ',m} bar{A}_{μ,m} ]

"""
function layer_pullback!(δ_prev::Array{<: AbstractGaugefields{NC,Dim},1},δ_current,layer::STOUT_Layer{Dim},Uprev,temps,tempf) where {NC,Dim}
    clear_U!(δ_prev)
    #δ_prev[μ](n) = δ_current[μ](n)*exp(Qμ[Uprev](n)) + F(δ_current,Uprev)
    #F(δ_current,Uprev) = sum_μ' sum_m Fm[μ](δ_current,Uprev)

    Cμs = similar(Uprev)
    construct_Cμ!(Cμs,layer,Uprev,temps)
    Qμs = similar(Uprev)
    construct_Qμs!(Qμs,Cμs,Uprev,temps)
    Λs = similar(Uprev)

    for μ=1:Dim
        construct_Λmatrix_forSTOUT!(Λs[μ],δ_current[μ],Qμs[μ],Uprev[μ])
        exptU!(δ_prev[μ],1,Qμs[μ],[temp1,temp2])
    end

    for μ=1:Dim
        for μd=1:Dim
            Λμd = Λs[μd]
            Uμd = Uprev[μd]

        end
    end

    temp1 = temps[1]
    temp2 = temps[2]

    for μ=1:Dim
        
    end

    F0 = tempf[1]
    temp1 = temps[1]
    Qμ = temps[2]


    
    Traceless_antihermitian_add!(Qμ,1,temp1)

    Λ = temps[2]

    #δ_next*exp(Q[Uprev]) 
    for μ=1:Dim
        #construct_Qμ!(Qμ,μ,Cμs,Uin,temps)
        construct_Λmatrix_forSTOUT!(Λ,δ_current,Qμs[μ],Uin[μ])


        mul!(temp1,Cμs[μ],Uprev[μ]') #Cμ*U^+
        clear_U!(F0)
        Traceless_antihermitian_add!(F0,1,temp1)
        exptU!(temp3,1,F0,[temp1,temp2])

        mul!(δ_next[μ],δ_prev[μ],temp3)   #δ_next*exp(Q[Uprev]) 
    end

end

"""
M = U δ star dexp(Q)/dQ
"""








function construct_Qμs!(Qμs,Cμs,Uin::Array{<: AbstractGaugefields{NC,Dim},1},temps) where {NC,Dim}
    temp1  = temps[1]
    for μ=1:Dim
        mul!(temp1,Cμs[μ],Uin[μ]') #Cμ*U^+
        clear_U!(Qμs)
        Traceless_antihermitian!(Qμs,temp1)
    end
end


function construct_Cμ!(Cμs,layer::STOUT_Layer{Dim},Uin::Array{<: AbstractGaugefields{NC,Dim},1},temps) where {NC,Dim}
    ρs = layer.ρs
    temp1  = temps[1]
    temp2  = temps[2]
    temp3  = temps[3]
    num = length(ρs)

    for μ=1:Dim
        clear_U!(Cμs[μ])
        for i=1:num
            loops = layer.dataset[i].Cμ[μ]
            evaluate_gaugelinks!(temp3,loops,Uin,[temp1,temp2])
            add_U!(Cμs[μ],ρs[i],temp3)
        end
    end
end

function apply_layer!(Uout::Array{<: AbstractGaugefields{NC,Dim},1},layer::STOUT_Layer{Dim},Uin,temps,tempf) where {NC,Dim}
    Cμs = similar(Uin)
    construct_Cμ!(Cμs,layer,Uin,temps)
    
    temp1  = temps[1]
    temp2  = temps[2]
    temp3  = temps[3]
    Qμ = tempf[1]

    for μ=1:Dim
        construct_Qμ!(Qμ,μ,Cμs,Uin,temps)
        # mul!(temp1,Cμs[μ],Uin[μ]') #Cμ*U^+
        #clear_U!(F0)
        #Traceless_antihermitian_add!(F0,1,temp1)
        exptU!(temp3,1,Qμ,[temp1,temp2])
        mul!(Uout[μ],temp3,Uin[μ])        
    end
    set_wing_U!(Uout)

    #=


    Cμ  = temps[4]
    temp1  = temps[1]
    temp2  = temps[2]
    temp3  = temps[3]

    F0 = tempf[1]
    ρs = layer.ρs

    num = length(ρs)

    for μ=1:Dim
        clear_U!(Cμ)
        for i=1:num
            loops = layer.dataset[i].Cμ[μ]
            evaluate_gaugelinks!(temp3,loops,Uin,[temp1,temp2])
            add_U!(Cμ,ρs[i],temp3)
        end
        mul!(temp1,Cμ,Uin[μ]') #Cμ*U^+
        clear_U!(F0)
        Traceless_antihermitian_add!(F0,1,temp1)
        
        exptU!(temp3,1,F0,[temp1,temp2])
        
        mul!(Uout[μ],temp3,Uin[μ])        
    end
    set_wing_U!(Uout)
    =#

end


function Stoutsmearing(loops_smearing,ρs)
    numloops = length(loops_smearing)
    staplesforsmear_set, tensor_derivative, staplesforsmear_dag_set, tensor_derivative_dag = make_loops_Stoutsmearing(loops_smearing,rand(numloops))
    return Stoutsmearing(
        staplesforsmear_set, 
        tensor_derivative,
        staplesforsmear_dag_set,
        tensor_derivative_dag,
        ρs
    )
end

function make_loops_Stoutsmearing(loops_smearing,ρs)
    num = length(loops_smearing)
    @assert num == length(ρs) "number of ρ elements in stout smearing scheme should be $num. Now $(length(ρs)). "
    staplesforsmear_set = Array{Wilson_loop_set,1}[]
    staplesforsmear_dag_set = Array{Wilson_loop_set,1}[]
    println("staple for stout smearing")

    tensor_derivative = Array{Tensor_derivative_set,1}(undef,num)
    tensor_derivative_dag = Array{Tensor_derivative_set,1}(undef,num)
    


    for i=1:num
        loop_smearing = loops_smearing[i]

        staplesforsmear = Array{Wilson_loop_set,1}(undef,4)
        staplesforsmear_dag = Array{Wilson_loop_set,1}(undef,4)
        


        staple = make_staples(loop_smearing)

        for μ=1:4
            staplesforsmear_dag[μ] = staple[μ]##make_plaq_staple(μ)
            staplesforsmear[μ] = staplesforsmear_dag[μ]'
            #println("$μ -direction")
            #display(staplesforsmear[μ])
            #staplesforsmear[μ] = make_plaq_staple(μ)
            #staplesforsmear_dag[μ] = staplesforsmear[μ]'
            #println("dagger: $μ -direction")
            #display(staplesforsmear_dag[μ])
        end
        tensor_derivative[i] = Tensor_derivative_set(staplesforsmear)
        tensor_derivative_dag[i] = Tensor_derivative_set(staplesforsmear_dag)

        push!(staplesforsmear_set,staplesforsmear )
        push!(staplesforsmear_dag_set,staplesforsmear_dag)

    end

    return staplesforsmear_set, 
        tensor_derivative, 
        staplesforsmear_dag_set, 
        tensor_derivative_dag 
        
end


function calc_coefficients_Q(Q)
    @assert size(Q) == (3,3)
    c0 = Q[1,1]*Q[2,2]*Q[3,3]+Q[1,2]*Q[2,3]*Q[3,1]+Q[1,3]*Q[2,1]*Q[3,2]-Q[1,3]*Q[2,2]*Q[3,1]-Q[1,2]*Q[2,1]*Q[3,3]-Q[1,1]*Q[2,3]*Q[3,2]
    #@time cdet = det(Q)
    ##println(c0,"\t",cdet)
    #exit() 
    
    c1 = 0.0
    for i=1:3
        for j=1:3
            c1 += Q[i,j]*Q[j,i]
        end
    end
    c1 /= 2
    c0max = 2*(c1/3)^(3/2)
    θ = acos(c0/c0max)
    u = sqrt(c1/3)*cos(θ/3)
    w = sqrt(c1)*sin(θ/3)
    ξ0 = sin(w)/w
    ξ1 = cos(w)/w^2 - sin(w)/w^3

    emiu = exp(-im*u)
    e2iu = exp(2*im*u)

    h0 = (u^2-w^2)*e2iu + emiu*(
        8u^2*cos(w)+2*im*u*(3u^2+w^2)* ξ0
    )
    h1 = 2u*e2iu-emiu*(
        2u*cos(w)-im*(3u^2-w^2)* ξ0
    )
    h2 = e2iu - emiu*(cos(w)+3*im*u*ξ0)

    denom = 9u^2-w^2
    
    f0 = h0/denom
    f1 = h1/denom
    f2 = h2/denom

    r10 = 2*(u+im*(u^2-w^2))*e2iu + 
            2*emiu*(
                4u*(2-im*u)*cos(w) + 
                im*(9u^2+w^2-im*u*(3u^2+w^2))*ξ0
            )
    r11 = 2*(1+2*im*u)*e2iu+ 
            emiu*(
                -2*(1-im*u)*cos(w)+
                im*(6u+im*(w^2-3u^2))*ξ0
            )
    r12 = 2*im*e2iu + im*emiu*(
        cos(w) -3*(1-im*u)*ξ0
    )
    r20 = -2*e2iu+2*im*u*emiu*(
        cos(w)+(1+4*im*u)*ξ0+3u^2*ξ1
    )
    r21 = -im*emiu*(
        cos(w)+(1+2*im*u)*ξ0 - 
        3*u^2*ξ1
    )
    r22 = emiu*(
        ξ0-3*im*u*ξ1
    )
    b10 = (
        2*u*r10+(3u^2-w^2)*r20 - 2*(15u^2+w^2)*f0
        )/(
            2*denom^2
        )
    
    b11 = (
            2*u*r11+(3u^2-w^2)*r21 - 2*(15u^2+w^2)*f1
            )/(
                2*denom^2
            )
    b12 = (
        2*u*r12+(3u^2-w^2)*r22 - 2*(15u^2+w^2)*f2
        )/(
            2*denom^2
        )
    b20 = (
        r10 - 3*u*r20 - 24*u*f0
        )/(2*denom^2)
    b21 = (
            r11 - 3*u*r21 - 24*u*f1
            )/(2*denom^2)
    b22 = (
        r12 - 3*u*r22 - 24*u*f2
        )/(2*denom^2)

    return f0,f1,f2,b10,b11,b12,b20,b21,b22
end

function calc_Bmatrix!(B,q,Q,NC)
    @assert NC == 2 "NC should be 2! now $NC"
    mul!(B,cos(q)/q -sin(q)/q^2,Q)
    for i=1:NC
        B[i,i] += -sin(q)
    end
    B .*= -1/2q
    #B[:,:] .= (cos(q)/q -sin(q)/q^2 )*Q

    #q = sqrt((-1/2)*tr(Q^2))
    #B = -(-sin(q)*I0_2 +(cos(q)/q -sin(q)/q^2 )*Q)*(1/2q)
end