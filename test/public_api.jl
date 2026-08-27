using Aqua
using LatticeQCD
using Test

@testset "Public API" begin
    exported_names = names(LatticeQCD; all=false, imported=false)
    @test !isempty(exported_names)
    for name in exported_names
        @test isdefined(LatticeQCD, name)
    end
end

@testset "Package quality" begin
    Aqua.test_all(
        LatticeQCD;
        ambiguities=false,
        persistent_tasks=false,
    )
    @test isempty(Test.detect_ambiguities(LatticeQCD; recursive=true))
end
