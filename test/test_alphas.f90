module test_cubic_alphas
    use yaeos__constants, only: pr
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use auxiliar_functions, only: allclose
    implicit none

    real(pr) :: absolute_tolerance = 1e-5_pr

contains
    subroutine collect_suite(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("AlphaSoave", test_alpha_soave), &
            new_unittest("AlphaRKPR", test_alpha_RKPR) &
            ]
    end subroutine collect_suite

    subroutine test_alpha_soave(error)
        use yaeos__constants, only: pr
        use yaeos__models_ar_cubic_alphas, only: AlphaSoave
        type(error_type), allocatable, intent(out) :: error

        type(AlphaSoave) :: alpha

        integer, parameter :: n = 2
        real(pr) :: Tr(n), k(n)
        real(pr) :: a(n), dadt(n), dadt2(n)

        real(pr) :: aval(n) = [1.1524213446238356, 1.0594365081389596]
        real(pr) :: davaldt(n) = [-0.33947331922020552, -0.14556349186104045]
        real(pr) :: davaldt2(n) = [0.47434164902525683, 0.15556349186104046]

        Tr = [0.4_pr, 0.5_pr]
        k = [0.2_pr, 0.1_pr]

        alpha%k = k
        call alpha%alpha(Tr, a, dadt, dadt2)

        call check(error, allclose(a, aval, absolute_tolerance))
        call check(error, allclose(dadt, davaldt, absolute_tolerance))
        call check(error, allclose(dadt2, davaldt2, absolute_tolerance))
    end subroutine test_alpha_soave

    subroutine test_alpha_RKPR(error)
        use yaeos__constants, only: pr
        use yaeos__models_ar_cubic_alphas, only: AlphaRKPR
        type(error_type), allocatable, intent(out) :: error

        type(AlphaRKPR) :: alpha

        integer, parameter :: n = 2
        real(pr) :: Tr(n), k(n)
        real(pr) :: a(n), dadt(n), dadt2(n)

        real(pr) :: aval(n) = [1.0456395525912732, 1.0183993761470242]
        real(pr) :: davaldt(n) = [-8.7136629382606107E-002, -4.0735975045880973E-002]
        real(pr) :: davaldt2(n) = [4.3568314691303053E-002, 1.7923829020187628E-002]

        Tr = [0.4_pr, 0.5_pr]
        k = [0.2_pr, 0.1_pr]

        alpha%k = k
        call alpha%alpha(Tr, a, dadt, dadt2)

        call check(error, allclose(a, aval, absolute_tolerance))
        call check(error, allclose(dadt, davaldt, absolute_tolerance))
        call check(error, allclose(dadt2, davaldt2, absolute_tolerance))
    end subroutine test_alpha_RKPR
end module test_cubic_alphas
