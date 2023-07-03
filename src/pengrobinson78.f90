module pengrobinson78
    !-| Peng Robinson 78 Equation of State.
    use yaeos_constants, only: pr, R
    use yaeos_autodiff
    use yaeos_ar_models, only: set_ar_function
    
    ! Use the generic Cubic EoS Ar function
    use yaeos_generic_cubic, only: set_functions, a_res

    ! Use the parameters defined in the Cubic EoS Module
    use yaeos_cubic_eos, only: set_parameters, &
                         a_classic, b_classic, c_classic, del1_classic, del2_classic
   
    ! Use the subroutines defined in the ClassicVdW module
    use yaeos_mixrule_classicvdw, only: setup_ClassicVdW, mix_ClassicVdW

    implicit none

    private
    public :: setup_pr78

contains
    subroutine setup_pr78(n, tc_in, pc_in, w_in, kij_in, lij_in)
        !-| Setup the enviroment to use the PengRobinson 78 Equation of State
        !   It uses the Cubic Van der Waals mixing rules
        integer :: n !! Number of components
        real(pr) :: tc_in(n)
        real(pr) :: pc_in(n)
        real(pr) :: w_in(n)
        real(pr) :: kij_in(n, n)
        real(pr) :: lij_in(n, n)

        real(pr) :: ac(n), b(n), c(n), k(n), del1(n), del2(n), tc(n), pc(n), w(n)

        del1 = 1._pr + sqrt(2._pr)
        del2 = 1._pr - sqrt(2._pr)

        tc = tc_in
        pc = pc_in
        w = w_in

        ! Setup the EoS parameters
        ac = 0.45723553_pr * R**2 * tc**2 / pc
        b = 0.07779607_pr * R * tc/pc
        where (w <= 0.491_pr)
           k = 0.37464 + 1.54226 * w - 0.26992 * w**2
        elsewhere
           k = 0.379642 + 1.48503 * w - 0.164423 * w**2 + 0.016666 * w**3
        end where

        c = 0

        ! Set the model parameters
        call set_parameters(ac, b, c, k, del1, del2, pc, tc, w)

        ! Setup the mixing rule parameters
        call setup_ClassicVdW(kij_in, lij_in)

        ! Setup the parameters functions in the generic cubicEoS module
        call set_functions(&
            a_classic, b_classic, c_classic, del1_classic, del2_classic, mix_ClassicVdW &
        )

        call set_ar_function(a_res)
    end subroutine

    subroutine setup_pr78_from_namelist(filepath)
        character(len=*), intent(in) :: filepath
        integer :: funit
        open(funit, file=filepath)
    end subroutine
end module
