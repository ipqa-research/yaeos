MODULE TapeRobinson
    use tapenade_pr
contains
    subroutine main()
        use yaeos_constants, only: pr
        integer, parameter :: n=2
        real(8) :: z(n), tc(n), pc(n), w(n), kij(n,n), lij(n,n)

        real(8) :: v, t

        real(8) :: ar, arv, arv2, art, art2, artv
        real(8) :: arn(n), arvn(n), artn(n), arn2(n,n)

        z = [0.3_pr, 0.7_pr]
        tc = [190._pr, 310._pr]
        pc = [14._pr, 30._pr]
        w = [0.001_pr, 0.03_pr]

        kij = reshape([0._pr, 0.1_pr, 0.1_pr, 0._pr], [n,n])
        lij = kij / 2._pr

        v = 1
        t = 150

        call SETUP_MODEL(tc, pc, w, kij, lij)
        call model%residual_helmholtz(&
            z, V, T, Ar=Ar, ArV=ArV, ArV2=ArV2, ArT=ArT, ArTV=ArTV, &
            ArT2=ArT2, Arn=Arn, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2 &
            )
        print *, "Ar: ", ar

        print *, "ArV: ", arV
        print *, "ArT: ", arT

        print *, "ArT2: ", arT2
        print *, "ArV2: ", ArV2

        print *, "ArTV: ", ArTV

        print *, "Arn: ", Arn

        print *, "ArVn: ", ArVn
        print *, "ArTn: ", ArTn

        print *, "Arn2: ", Arn2
    end subroutine main
end module TapeRobinson
