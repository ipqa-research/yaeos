module test_cubic_mixrules
    use yaeos_constants, only: pr
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use auxiliar_functions, only: allclose
    implicit none

    real(pr) :: absolute_tolerance = 1e-5_pr

contains
    subroutine collect_suite(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("QMR_RKPR", test_QMR_RKPR) &
            ]
    end subroutine collect_suite

    subroutine test_QMR_RKPR(error)
        use yaeos_constants, only: pr
        use yaeos_models_ar_cubic_quadratic_mixing, only: QMR_RKPR
        type(error_type), allocatable, intent(out) :: error

        type(QMR_RKPR) :: mixrule

        integer, parameter :: n = 3
        real(pr) :: del1(3) = [0.2, 0.5, 0.6]
        real(pr) :: z(3), D1, dD1i(n), dD1ij(n, n)
        real(pr) :: d1_val, dd1i_val(n), dd1ij_val(n, n)

        z = [0.3, 0.5, 0.2]

        call mixrule%D1mix(z, del1, D1, dD1i, dD1ij)
         d1_val = 0.43000000342726713     
        dD1i_val = [-0.22999999701976787, 6.9999995529651651E-002,  0.17000001788139310]
        dD1ij_val = reshape([0.45999998718500179, 0.15999999910593043, 5.9999978244305405E-002, &
                             0.15999999910593043, -0.13999998897314089, -0.24000000983476594, &
                             5.9999978244305405E-002, -0.24000000983476594,      -0.34000003069639095], [n,n])

        call check(error, allclose([D1], [D1_val], absolute_tolerance))
        call check(error, allclose([dD1i], [dD1i_val], absolute_tolerance))
        call check(error, allclose([dD1ij], [dD1ij_val], absolute_tolerance))
        

    end subroutine test_QMR_RKPR

end module test_cubic_mixrules
