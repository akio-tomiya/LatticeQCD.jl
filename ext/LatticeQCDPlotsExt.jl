module LatticeQCDPlotsExt

using LatticeQCD
using Plots
import LatticeQCD: demo

function demo_plot(plaquettes, polyakovs)
    trajectories = eachindex(plaquettes)
    magnitudes = abs.(polyakovs)
    bins = max(1, round(Int, 1.5 * log(length(plaquettes)) + 1))
    return plot(
        histogram(plaquettes; bins, xlabel="Plaquette", title="SU(3), quenched"),
        plot(trajectories, plaquettes; xlabel="MC time", ylabel="Plaquette",
             title="L=4^4, Heatbath"),
        histogram(magnitudes; bins, xlabel="|Polyakov loop|"),
        plot(trajectories, magnitudes; xlabel="MC time", ylabel="|Polyakov loop|"),
        scatter(real.(polyakovs), imag.(polyakovs); xlabel="Re", ylabel="Im",
                title="Polyakov loop"),
        plot(trajectories, angle.(polyakovs); xlabel="MC time",
             ylabel="Arg(Polyakov loop)"),
        layout=(3, 2), legend=false,
    )
end

function demo(; steps::Integer=10_000, display_every::Integer=1, showplot::Bool=true)
    steps > 0 || throw(ArgumentError("steps must be positive"))
    display_every > 0 || throw(ArgumentError("display_every must be positive"))
    config = LQCDConfig(
        LatticeConfig((4, 4, 4, 4)),
        GaugeConfig(3, 1, "hot", nothing, 111),
        GaugeActionConfig(
            GaugeActionTermConfig(:gauge_plaquette, "plaquette", 6.758506610326683),
        ),
        HeatbathConfig(true, 100_000, 0, RandomStreamConfig(111, :heatbath)),
    )
    program = MeasurementProgram(MeasurementPlan(
        ScheduledMeasurementConfig(PlaquetteObservableConfig(), PeriodicSchedule(1)),
        ScheduledMeasurementConfig(PolyakovLoopObservableConfig(), PeriodicSchedule(1)),
    ))
    schedule = SimulationSchedule(0, steps; measurements=program)
    session = build_simulation(SimulationSpec(config, schedule, OutputConfig()))
    plaquettes = Float64[]
    polyakovs = ComplexF64[]
    figure = nothing
    for trajectory in 1:steps
        records = step!(session).measurements
        push!(plaquettes, records[1].value)
        push!(polyakovs, records[2].value)
        if trajectory == steps || (showplot && mod(trajectory, display_every) == 0)
            figure = demo_plot(plaquettes, polyakovs)
            showplot && display(figure)
        end
    end
    return figure
end

end
