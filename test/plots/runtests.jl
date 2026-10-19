using Test
import Logging

ENV["GKSwstype"] = "100"

@testset "Plots loading order" begin
    if only(ARGS) == "lattice-first"
        @eval using LatticeQCD
        @test Base.get_extension(LatticeQCD, :LatticeQCDPlotsExt) === nothing
        @test_logs min_level=Logging.Warn @eval using Plots
    elseif only(ARGS) == "plots-first"
        @eval using Plots
        @test_logs min_level=Logging.Warn @eval using LatticeQCD
    else
        error("expected lattice-first or plots-first")
    end
    @test Base.get_extension(LatticeQCD, :LatticeQCDPlotsExt) !== nothing
    @test applicable(demo)
end

mutable struct DemoTestDisplay <: AbstractDisplay
    frames::Int
end
Base.display(display::DemoTestDisplay, ::Plots.Plot) = (display.frames += 1; nothing)

@testset "Heatbath demo draws and returns a plot" begin
    @test_throws ArgumentError demo(steps=0)
    @test_throws ArgumentError demo(display_every=0)
    display = DemoTestDisplay(0)
    pushdisplay(display)
    try
        mktempdir() do directory
            cd(directory) do
                figure = demo(steps=3, showplot=false)
                @test figure isa Plots.Plot
                @test length(figure.subplots) == 6
                @test display.frames == 0
                @test isempty(readdir())
                plaquettes = figure.subplots[2].series_list[1][:y]
                @test length(plaquettes) == 3
                @test all(isfinite, plaquettes)
                savefig(figure, "demo.png")
                @test filesize("demo.png") > 100
                @test read("demo.png")[1:8] == UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]
                refreshed = demo(steps=3, display_every=2)
                @test display.frames == 2 # interval and final frame
                @test refreshed.subplots[2].series_list[1][:y] == plaquettes
            end
        end
    finally
        popdisplay(display)
    end
end
