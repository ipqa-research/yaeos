module tapenade_model_template
    use yaeos_tapenade_ar_api, only: ArModelTapenade
    implicit none

    real(8), parameter :: r=0.08314472 ! < keep the R constant used in the module
    type(ArModelTapenade) :: model

    <
        add your model parameters here
    >

contains
    subroutine setup_model(<your set of parameters>)
        <
            set of paremeters
        >

        <
            Setup the module parameters from the input 
        >

        ! This must be kept to setup the tapemodel object
        model%ar => ar
        model%ar_d => ar_d
        model%ar_b => ar_b
        model%ar_d_b => ar_d_b
        model%ar_d_d => ar_d_d
    end subroutine

    subroutine ar(n, v, t, arval)
        real(8), intent(in) :: n(:), v, t
        real(8), intent(out) :: arval
        ! ^ These names should not be modified

        <
            any other variables definitions
        >

        <
            Calculate arval
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
