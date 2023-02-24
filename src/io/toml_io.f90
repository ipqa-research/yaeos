module toml_io
    !! This module handles the IO of the input files
    !!
    !! The inputs files are based on a TOML structure, where there are two main
    !! sets of tables:
    !! - [system table]
    !! - [component tables]
    !!
    !! The system table handles all the main configurations, like thermodynamic
    !! model to use, names of components, concentrations, etc
    !!
    !! The component tables are a set of subtables for each component, where
    !! each table correspond to a specific model, with the exception of the
    !! "<component>.pure" table, which holds the critical constants of the
    !! component
    !!
    !! Example simple configuration file
    !! ```toml
    !! [system]
    !! model = "SoaveRedlichKwong"
    !! use_parameters = false  # Use EOS parameters instead of critical constants
    !! components = [ "methane", "ethane" ]
    !! composition = [ 0.3, 0.6 ]
    !! mixing_rule = "ClassicVdW"
    !! 
    !! [methane.pure]
    !! name = "C1"
    !! Tc = 190.6000
    !! Pc = 46.0000
    !! Vc = 0.008000
    !! w = 0.11484
    !! 
    !! [methane.SoaveRedlichKwong]
    !! a = 2.3339
    !! b = 0.029848
    !! k = 0.492581
    !! 
    !! # ethane tables can be located inside an external ethane.toml file
    !! ```
    use thermo_io
    use toml_cubic


    implicit none


    private

    public :: read_system, size!, setup_from_toml

end module toml_io
