program test_alphas
    use yaeos__constants, only: pr
    use yaeos__models_ar_cubic_alphas, only: AlphaSoave, AlphaRKPR
    use auxiliar_functions, only: allclose
    use testing_aux, only: assert, test_title
    implicit none

    real(pr) :: absolute_tolerance = 1e-5_pr

    write(*, *) test_title("CUBIC ALPHA FUNCTIONS TESTS")

    call test_alpha_soave()
    call test_alpha_RKPR()

    write(*, *) " "

contains
    subroutine test_alpha_soave()
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

        call assert(allclose(a, aval, absolute_tolerance), "AlphaSoave: a values")
        call assert(allclose(dadt, davaldt, absolute_tolerance), "AlphaSoave: dadt values")
        call assert(allclose(dadt2, davaldt2, absolute_tolerance), "AlphaSoave: dadt2 values")
    end subroutine test_alpha_soave

    subroutine test_alpha_RKPR()
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

        call assert(allclose(a, aval, absolute_tolerance), "AlphaRKPR: a values")
        call assert(allclose(dadt, davaldt, absolute_tolerance), "AlphaRKPR: dadt values")
        call assert(allclose(dadt2, davaldt2, absolute_tolerance), "AlphaRKPR: dadt2 values")
    end subroutine test_alpha_RKPR
end program test_alphas
