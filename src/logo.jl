function print_wizard_logo(outs)
    blue    = "\033[34m"
    red     = "\033[31m"
    green   = "\033[32m"
    magenta = "\033[35m"
    normal  = "\033[0m\033[0m"

    logo = raw"""
--------------------------------------------------------------------------------  
run_wizard

LatticeQCD.jl                                               
    """



    logo = replace(logo, "Q" => "$(red)Q$(normal)")
    logo = replace(logo, "C" => "$(blue)C$(normal)")
    logo = replace(logo, "D" => "$(green)D$(normal)")
    println(outs, logo)
    println(outs,
        "Welcome to the LatticeQCD wizard.  ",
        "We'll get you set up in no time."
    )
    println("--------------------------------------------------------------------------------")
end
print_wizard_logo(stdout)