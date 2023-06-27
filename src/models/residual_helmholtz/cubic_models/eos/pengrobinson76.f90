module pengrobinson76
    use constants, only: pr, R
    use hyperdual_mod
    use generic_cubic, only: set_parameters_functions, set_mixrule
    use cubic_eos, only: ac, b, c, del1, del2, k, &
                         tc, pc, w, &
                         alloc, &
                         a_classic, b_classic, c_classic, del1_classic, del2_classic
    use mixrule_classicvdw, only: setup_ClassicVdW, mix_ClassicVdW

    implicit none

    private
    public :: setup_pr76

contains

    subroutine setup_pr76(n, tc_in, pc_in, w_in, kij_in, lij_in)
        ! use stdlib_optval, only: optval
        integer :: n
        real(pr) :: tc_in(n), pc_in(n), w_in(n)
        real(pr) :: kij_in(n, n), lij_in(n, n)

        del1 = 1._pr + sqrt(2.0_pr)
        del2 = 1._pr - sqrt(2.0_pr)

        tc = tc_in
        pc = pc_in
        w = w_in

        ac = 0.45723553_pr * R**2 * tc**2 / pc
        b = 0.07779607_pr * R * tc/pc
        k = 0.37464_pr + 1.54226_pr * w - 0.26993_pr * w**2
        c = 0

        call setup_ClassicVdW(kij_in, lij_in)
        call set_mixrule(mix_ClassicVdW)
        call set_parameters_functions(&
            a_classic, b_classic, c_classic, del1_classic, del2_classic &
        )
    end subroutine

end module
