["HMC related"]
"Δτ" = 0.05
MDsteps = 20

["Physical setting(fermions)"]
quench = false
Dirac_operator = "Wilson"

["Physical setting"]
L = [4, 4, 4, 4]
update_method = "HMC"
"β" = 5.7

["System Control"]
logfile = "HMC_L04040404_beta5.7_Wilson_kappa0.141139.txt"
log_dir = "./logs"

["Measurement set"]
measurement_dir = "HMC_L04040404_beta5.7_Wilson_kappa0.141139"
measurement_basedir = "./measurements"

    ["Measurement set".measurement_methods.Pion_correlator]
    measure_every = 1
    methodname = "Pion_correlator"

        ["Measurement set".measurement_methods.Pion_correlator.fermion_parameters]
        Dirac_operator = "Wilson"

    ["Measurement set".measurement_methods.Polyakov_loop]
    measure_every = 1
    methodname = "Polyakov_loop"

    ["Measurement set".measurement_methods.Plaquette]
    measure_every = 1
    methodname = "Plaquette"

    ["Measurement set".measurement_methods.Energy_density]
    measure_every = 1
    methodname = "Energy_density"
