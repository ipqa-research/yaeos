module tape_nrtl
    use yaeos_constants, only: pr
    use yaeos_models_ge, only: GeModel
    use yaeos_models_ge_NRTL, only: NRTL
contains
    subroutine main
        class(GeModel), allocatable :: model
        integer, parameter :: n=2
        real(pr) :: a(n, n), b(n, n), c(n, n)

        real(pr) :: ge, GeT, GeT2, Gen(n), Gen2(n, n), GeTn(n)
        real(pr) :: lngamma(n)

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

        call model%excess_gibbs(&
            moles, t, &
            Ge=Ge, GeT=GeT, Get2=GeT2, Gen=Gen, GeTn=GeTn, Gen2=Gen2 &
        )

        print *, "Ge: ", Ge
        print *, "GeT: ", GeT
        print *, "GeT2: ", GeT2
        print *, "Gen: ", Gen
        print *, "GeTn: ",  GeTn
        print *, "Gen2: ",  Gen2

        call model%ln_activity_coefficient(moles, T, lngamma)
        print *, exp(lngamma)

    end subroutine
end module