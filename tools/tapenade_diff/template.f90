module tapenade_model_template
    use yaeos__tapenade_ar_api, only: ArModelTapenade
    use yaeos__tapenade_interfaces
    use yaeos, only: R !! Ideal gas constants used on yaeos
    implicit none

    ! The comments on the type definition are what must be fixed by hand after
    ! differentiation, see the fixed result at tapeout/pr_diff.f90 for an 
    ! already fixed example

    !type, extends(ArModelTapenade) :: ModelName
    type :: ModelName
        <
            add your model parameters here
            ! For example
            ! real(8), allocatable :: critical_temperature(:)
        >
    ! contains
    !  procedure :: ar
    !  procedure :: ar_d
    !  procedure :: ar_b
    !  procedure :: ar_d_b
    !  procedure :: ar_d_d
    !  procedure :: v0
    end type ModelName


contains
    type(ModelName) function setup_model(<your set of parameters>)
        <
            set of paremeters
            ! for example:
            ! real(8), intent(in) :: tc(:)
        >

        <
            Setup the module parameters from the input 
            ! For example:
            ! setup_model%critical_temperature = tc
        >
    end function

    subroutine ar(model, n, v, t, arval)
        type(ModelName), intent(in) :: model
        real(8), intent(in) :: n(:), v, t
        real(8), intent(out) :: arval
        ! ^ These names should not be modified
        <
            any other variables definitions
            ! All the models parameters must be defined as internal variables
            ! even if they are already defined in the derived type, this
            ! assures that the differentiator gets their dimension related
            ! to `n`
            ! real(8) :: tc(size(n))
            ! later on we extract the value from the type:
            ! tc = model%critical_temperature
        >

        <
            Calculate arval
            ! Here goes Ar(n, V, T)
        >
    end subroutine

    pure function volume_initalizer(n, p, t) result(v0)
         real(8), intent(in) :: n(:)
         real(8), intent(in) :: p
         real(8), intent(in) :: t
         real(8) :: v0
        ! ^ These names should not be modified
         v0 = <how you initialize liquid volume for a pressure newton solver>
    end function
end module
