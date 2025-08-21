program test_legacy
    use yaeos__constants, only: pr, r
    use legacy_ar_models, only: setup, fact => PR76_factory, ar_srkpr, ArVnder, &
        kij, lij, ac, b, wmod => w, k
    use legacy_thermo_properties, only: termo
    use testing_aux, only: assert, test_title
    implicit none

    write(*, *) test_title("LEGACY MODELS TESTS")

    call test_legacy_ar()
    call test_legacy_termo()

    write(*, *) " "

contains
    subroutine test_legacy_ar()
        integer, parameter :: n = 2
        real(pr) :: tc(n), pc(n), w(n), z(n)
        real(pr) :: t, p, v
        real(pr) :: ar, arv, artv, arv2, arn(n), arvn(n), artn(n), arn2(n, n), et, st

        z = [0.3, 0.7]
        tc = [190, 310]
        pc = [14, 30]
        w = [0.001, 0.03]

        call setup(n, 1, 0, 1)
        call fact(z, tc_in=tc, pc_in=pc, w_in=w)

        kij = 0
        lij = 0

        v = 1.0
        t = 150

        call ArVnder(n, 2, 1, z, v, t, ar, arv, artv, arv2, arn, arvn, artn, arn2)

        call assert(maxval(abs(arn - [-15.996794004778227, -20.3860725744])) < 1e-5, "legacy_AR: arn")
        call assert(maxval(abs(arvn - [14.047950956798386, 18.367895313850489])) < 1e-5, "legacy_AR: arvn")
        call assert(maxval(abs(artn - [5.0336216865122976E-002, 5.1916228935932923E-002])) < 1e-5, "legacy_AR: artn")
        call assert(maxval(abs(arn2(1, :) - [-11.441418892682893, -15.165036190546420])) < 1e-5, "legacy_AR: arn2(1,:)")
        call assert(maxval(abs(arn2(2, :) - [-15.165036190546420, -19.740549301758474])) < 1e-5, "legacy_AR: arn2(2,:)")
    end subroutine test_legacy_ar

    subroutine test_legacy_termo()
        integer, parameter :: n = 2
        real(pr) :: tc(n), pc(n), w(n), z(n)
        real(pr) :: lnphi_p(n), lnphi_nn(n, n)
        real(pr) :: t, p, v

        z = [0.3, 0.7]
        tc = [190, 310]
        pc = [14, 30]
        w = [0.001, 0.03]

        call setup(n, 1, 0, 1)
        call fact(z, tc_in=tc, pc_in=pc, w_in=w)

        kij = 0
        lij = 0

        p = 1.0
        t = 150

        call termo(n, 1, 4, t, p, z, v, philog=lnphi_p, fugn=lnphi_nn)
        call assert(maxval(abs(lnphi_p - [1.7282684886888058, -2.3291258019214967])) < 1e-5, "legacy_TERMO: lnphi_p")
    end subroutine test_legacy_termo
end program test_legacy
