using LatticeQCD, Test

@testset "Plots remains optional" begin
    @test Base.find_package("Plots") === nothing
    @test all(id.name != "Plots" for id in keys(Base.loaded_modules))
    @test Base.get_extension(LatticeQCD, :LatticeQCDPlotsExt) === nothing
    @test :demo in names(LatticeQCD)
    @test !applicable(demo)
end
