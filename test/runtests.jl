using LatticeQCD
using Test


@testset "LatticeQCD.jl" begin
    eps = 1e-5
    fout = open("Testvalues.txt","w");println(fout,"Test values")
    @time begin # time
    @testset "quenched HMC" begin

        @testset "SU(2)" begin
            @time plaq = run_LQCD("test02.jl")
            #plaq_comparison = 0.5575491312570713
            plaq_comparison = 0.48785061558469406
            @test abs(plaq - plaq_comparison)/plaq_comparison < eps
            println(fout,"qhmc SU(2), $plaq")
        end

        @testset "SU(3)" begin
            @time plaq = run_LQCD("test01.jl")
            #plaq_comparison = 0.6190393357419764
            plaq_comparison = 0.5753703885492326
            @test abs(plaq - plaq_comparison)/plaq_comparison < eps
            println(fout,"qhmc SU(3), $plaq")
        end


        @testset "SU(4)" begin
            @time plaq = run_LQCD("test03.jl")
            #plaq_comparison = 0.4966683811089479
            plaq_comparison = 0.3915274700011152
            @test abs(plaq - plaq_comparison)/plaq_comparison < eps
            println(fout,"qhmc SU(4), $plaq")
        end
    end
    
    
    
    @testset "Heatbath" begin
        @testset "SU(2)" begin
            @time plaq = run_LQCD("test02-hb.jl")
            #plaq_comparison = 0.5287735727118359
            plaq_comparison = 0.4855523668804699
            @test abs(plaq - plaq_comparison)/plaq_comparison < eps
            println(fout,"hb SU(2), $plaq")
        end

        @testset "SU(3)" begin
            @time plaq = run_LQCD("test01-hb.jl")
            #plaq_comparison = 0.5821680570717788
            plaq_comparison = 0.5502269475635925
            @test abs(plaq - plaq_comparison)/plaq_comparison < eps
            println(fout,"hb SU(3), $plaq")
        end


        @testset "SU(4)" begin
            @time plaq = run_LQCD("test03-hb.jl")
            #plaq_comparison = 0.5467724338528576
            plaq_comparison = 0.4425954597477664
            @test abs(plaq - plaq_comparison)/plaq_comparison < eps
            println(fout,"hb SU(4), $plaq")
        end
    end



    
        
    @testset "HMC" begin

        @testset "Wilson SU(3) with SextonWeingargten" begin
            @time plaq = run_LQCD("test_wilson.jl")
            #plaq_comparison = 0.3449688128155864
            #plaq_comparison = 0.4867607208994073
            plaq_comparison = 0.4976172009730353
            @test abs(plaq - plaq_comparison)/plaq_comparison < eps
            println(fout,"Wilson SU(3) with SextonWeingargten, $plaq")
        end

        @testset "Staggered SU(3) with 4 tastes" begin
            @time plaq = run_LQCD("test_staggered.jl")
            #plaq_comparison = 0.25455870400018477 #0.00624053999484795
            #plaq_comparison = 0.4987738124715037
            plaq_comparison = 0.4713469212809392
            @test abs(plaq - plaq_comparison)/plaq_comparison < eps
            println(fout,"Staggered SU(3) with 4 tastes, $plaq")
        end

        @testset "Staggered SU(3) with 2 tastes" begin
            @time plaq = run_LQCD("test_Nf2.jl")
            #plaq_comparison = 0.5630198767336069
            plaq_comparison = 0.5837848292310798
            @test abs(plaq - plaq_comparison)/plaq_comparison < eps
            println(fout,"Staggered SU(3) with 2 tastes, $plaq")
        end

        @testset "Staggered SU(3) with 3 tastes" begin
            @time plaq = run_LQCD("test_Nf3.jl")
            #plaq_comparison = 0.565176584402352
            plaq_comparison = 0.5864438294310259
            @test abs(plaq - plaq_comparison)/plaq_comparison < eps
            println(fout,"Staggered SU(3) with 3 tastes, $plaq")
        end

        @testset "Domain-wall SU(3)" begin
            # cold start
            @time plaq = run_LQCD("test_domainwallhmc.toml")
            plaq_comparison = 0.649479122018118
            @test abs(plaq - plaq_comparison)/plaq_comparison < eps
            println(fout,"Domain-wall SU(3), $plaq")
        end
    end
    end # time
    close(fout)

    #=       
    @testset "SLMC" begin

        @testset "Staggered SU(3)" begin
            @time plaq = run_LQCD("test06_slmc_ks.jl")
            plaq_comparison = 0.5279623287035208
            @test abs(plaq - plaq_comparison)/plaq_comparison < eps
        end
    end
        @testset "WilsonClover SU(3) with SextonWeingargten" begin
            @time plaq = run_LQCD("test_wilsonclover.jl")
            plaq_comparison = 0.3449688128155864
            #plaq_comparison = 0.00624053999484795
            @test abs(plaq - plaq_comparison)/plaq_comparison < eps
        end
    =#

    
    

end

