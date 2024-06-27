module test_legacy
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use auxiliar_functions
    implicit none

contains
    subroutine collect_suite(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("legacy_AR", test_legacy_ar), &
            new_unittest("legacy_TERMO", test_legacy_termo) &
            ]
    end subroutine collect_suite

    subroutine test_legacy_ar(error)
        use yaeos__constants, only: pr, r
        use legacy_ar_models, only: setup, fact => PR76_factory, ar_srkpr, ArVnder, &
            kij, lij, ac, b, wmod => w, k

        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 2
        real(pr) :: tc(n), pc(n), w(n), z(n)
        real(pr) :: t, p, v
        real(pr) :: ar, arv, artv, arv2, arn(n), arvn(n), artn(n), arn2(n, n), et, st

        integer :: i

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

        ! call check(error, abs(ar - -9.8908394285700) < 1e-5)
        ! call check(error, abs(arv - 9.1784903035682) < 1e-5)
        ! call check(error, abs(artv - -2.52238310493E-002) < 1e-5)
        ! call check(error, abs(arv2 - -17.0719465497) < 1e-5)
        call check(error, maxval(abs(arn - [-15.996794004778227, -20.3860725744])) < 1e-5)
        call check(error, maxval(abs(arvn - [14.047950956798386, 18.367895313850489])) < 1e-5)
        call check(error, maxval(abs(artn - [5.0336216865122976E-002, 5.1916228935932923E-002])) < 1e-5)
        call check(error, maxval(abs(arn2(1, :) - [-11.441418892682893, -15.165036190546420])) < 1e-5)
        call check(error, maxval(abs(arn2(2, :) - [-15.165036190546420, -19.740549301758474])) < 1e-5)
    end subroutine test_legacy_ar

    subroutine test_legacy_termo(error)
        use yaeos__constants, only: pr, r
        use legacy_ar_models, only: setup, fact => PR76_factory, &
            kij, lij, ac, b, wmod => w, k
        use legacy_thermo_properties, only: termo

        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 2
        real(pr) :: tc(n), pc(n), w(n), z(n)
        real(pr) :: lnphi_p(n), lnphi_nn(n, n)
        real(pr) :: t, p, v

        integer :: i

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
        call check(error, maxval(abs(lnphi_p - [1.7282684886888058, -2.3291258019214967])) < 1e-5)
    end subroutine test_legacy_termo
end module test_legacy
