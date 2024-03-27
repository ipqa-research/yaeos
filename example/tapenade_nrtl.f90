module tape_nrtl
    use yaeos_constants, only: pr
    use yaeos_models_ge, only: GeModel
    use yaeos_models_ge_NRTL, only: NRTL
contains
    subroutine main
        class(GeModel), allocatable :: model
        integer, parameter :: n=2
        real(pr) :: a(n, n), b(n, n), c(n, n)

        real(pr) :: ge, dgedt, dgedn(n)

        real(pr) :: T = 150
        real(pr) :: moles(n) = [0.3, 0.7]

        a = 0; b = 0; c = 0

        a(1, 2) = 3.458
        a(2, 1) = -0.801

        b(1, 2) = -586.1
        b(2, 1) = 246.2

        c(1, 2) = 0.3
        c(2, 1) = 0.3

        model = NRTL(a, b, c)

        call model%excess_gibbs(moles, t, ge=ge, GeT=dgedt, Gen=dgedn)
        print *, ge, dgedt
        call model%ln_activity_coefficient(moles, t, dgedn)
        print *, dgedn
        print *, exp(dgedn)

    end subroutine
end module