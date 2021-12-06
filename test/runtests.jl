using LatticeQCD
using Test


@testset "LatticeQCD.jl" begin
    eps = 1e-1
    
    
    @testset "quenched HMC" begin

        @testset "SU(2)" begin
            @time plaq = run_LQCD("test02.jl")
            plaq_comparison = 0.5575491312570713
            @test abs(plaq - plaq_comparison)/plaq_comparison < eps
        end

        @testset "SU(3)" begin
            @time plaq = run_LQCD("test01.jl")
            plaq_comparison = 0.6190393357419764
            @test abs(plaq - plaq_comparison)/plaq_comparison < eps
        end


        @testset "SU(4)" begin
            @time plaq = run_LQCD("test03.jl")
            plaq_comparison = 0.4966683811089479
            @test abs(plaq - plaq_comparison)/plaq_comparison < eps
        end
    end
    
    
    
    @testset "Heatbath" begin
        @testset "SU(2)" begin
            @time plaq = run_LQCD("test02-hb.jl")
            plaq_comparison = 0.5287735727118359
            @test abs(plaq - plaq_comparison)/plaq_comparison < eps
        end

        @testset "SU(3)" begin
            @time plaq = run_LQCD("test01-hb.jl")
            plaq_comparison = 0.5821680570717788
            @test abs(plaq - plaq_comparison)/plaq_comparison < eps
        end


        @testset "SU(4)" begin
            @time plaq = run_LQCD("test03-hb.jl")
            plaq_comparison = 0.5467724338528576
            @test abs(plaq - plaq_comparison)/plaq_comparison < eps
        end
    end



    
        
    @testset "HMC" begin

        @testset "WilsonClover SU(3) with SextonWeingargten" begin
            @time plaq = run_LQCD("test_wilsonclover.jl")
            plaq_comparison = 0.3449688128155864
            #plaq_comparison = 0.00624053999484795
            @test abs(plaq - plaq_comparison)/plaq_comparison < eps
        end

        @testset "Staggered SU(3) with 4 tastes" begin
            @time plaq = run_LQCD("test_staggered.jl")
            plaq_comparison = 0.25455870400018477 #0.00624053999484795
            @test abs(plaq - plaq_comparison)/plaq_comparison < eps
        end

        @testset "Staggered SU(3) with 2 tastes" begin
            @time plaq = run_LQCD("test_Nf2.jl")
            plaq_comparison = 0.5630198767336069
            @test abs(plaq - plaq_comparison)/plaq_comparison < eps
        end

        @testset "Staggered SU(3) with 3 tastes" begin
            @time plaq = run_LQCD("test_Nf3.jl")
            plaq_comparison = 0.565176584402352
            @test abs(plaq - plaq_comparison)/plaq_comparison < eps
        end
    end

            
    @testset "SLMC" begin

        @testset "Staggered SU(3)" begin
            @time plaq = run_LQCD("test06_slmc_ks.jl")
            plaq_comparison = 0.5279623287035208
            @test abs(plaq - plaq_comparison)/plaq_comparison < eps
        end
    end

    
    

end

