! In `yaeos` there are routines for calculation of saturation points. To obtain
! either bubble, dew or liquid-liquid points
!
! - saturation_pressure
! - saturation_temperature
! 
! All the outputs of this functions use the `EquilibriumState` type
program saturation
    use yaeos
    use yaeos__example_tools, only: methane_butane_pr76

    class(ArModel), allocatable :: model
    type(EquilibriumState) :: sat_point

    real(pr) :: n(2), T

    ! This function returns a predefined model to simplify the step. How it 
    ! works can be seen in the `example_tools.f90` file.
    model = methane_butane_pr76()

    n = [2.5, 6.7]

    ! ==========================================================================
    ! Calculate the bubble pressure of this system
    ! --------------------------------------------------------------------------
    write(*, *) "Bubble pressure:"
    T = 150
    sat_point = saturation_pressure(model, n, T=T, kind="bubble")
    write (*, *) "kind, T, P: ", sat_point%kind, sat_point%T, sat_point%P
    write (*, *) "x: ", sat_point%x
    write (*, *) "y: ", sat_point%y
    
    ! ==========================================================================
    ! Calculate the bubble temperature of this system
    ! --------------------------------------------------------------------------
    write(*, *) ""
    write (*, *) "Bubble temperature:"
    sat_point = saturation_temperature(model, n, P=15._pr, kind="bubble")
    write (*, *) "kind, T, P: ", sat_point%kind, sat_point%T, sat_point%P
    write (*, *) "x: ", sat_point%x
    write (*, *) "y: ", sat_point%y

    ! ==========================================================================
    ! Calculate the dew pressure of this system
    ! --------------------------------------------------------------------------
    write(*, *) ""
    write(*, *) "Dew pressure:"
    sat_point = saturation_pressure(model, n, T=150._pr, kind="dew")
    write (*, *) "kind, T, P: ", sat_point%kind, sat_point%T, sat_point%P
    write (*, *) "x: ", sat_point%x
    write (*, *) "y: ", sat_point%y
    
    ! ==========================================================================
    ! Calculate the dew temperature of this system
    ! --------------------------------------------------------------------------
    write(*, *) ""
    write (*, *) "Dew temperature:"
    sat_point = saturation_temperature(model, n, P=15._pr, kind="dew", t0=330._pr)
    write (*, *) "kind, T, P: ", sat_point%kind, sat_point%T, sat_point%P
    write (*, *) "x: ", sat_point%x
    write (*, *) "y: ", sat_point%y
end program