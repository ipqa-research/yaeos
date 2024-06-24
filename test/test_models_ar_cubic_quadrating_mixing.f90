module test_cubic_mixrules
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
                  new_unittest("QMR_RKPR", test_QMR_RKPR), &
                  new_unittest("MHV_SRK", test_MHV_SRK) &
                  ]
   end subroutine collect_suite

   subroutine test_QMR_RKPR(error)
      use yaeos__constants, only: pr
      use yaeos__models_ar_cubic_quadratic_mixing, only: QMR_RKPR
      type(error_type), allocatable, intent(out) :: error

      type(QMR_RKPR) :: mixrule

      integer, parameter :: n = 3
      real(pr) :: del1(3) = [0.2, 0.5, 0.6]
      real(pr) :: z(3), D1, dD1i(n), dD1ij(n, n)
      real(pr) :: d1_val, dd1i_val(n), dd1ij_val(n, n)

      z = [0.3, 0.5, 0.2]

      call mixrule%D1mix(z, del1, D1, dD1i, dD1ij)
      d1_val = 0.43000000342726713
      dD1i_val = [-0.22999999701976787, 6.9999995529651651E-002, 0.17000001788139310]
      dD1ij_val = reshape([0.45999998718500179, 0.15999999910593043, 5.9999978244305405E-002, &
                           0.15999999910593043, -0.13999998897314089, -0.24000000983476594, &
                           5.9999978244305405E-002, -0.24000000983476594, -0.34000003069639095], [n, n])

      call check(error, allclose([D1], [D1_val], absolute_tolerance))
      call check(error, allclose([dD1i], [dD1i_val], absolute_tolerance))
      call check(error, allclose([dD1ij], [dD1ij_val], absolute_tolerance))

   end subroutine test_QMR_RKPR

   subroutine test_MHV_SRK(error)
      use yaeos__constants, only: pr
      use yaeos__models, only: CubicEoS, NRTL, MHV
      use fixtures_models, only: binary_NRTL_SRK
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 2
      real(pr) :: test_D = 30.797085158580309
      real(pr) :: test_dDdT = -8.4060598003984383E-002
      real(pr) :: test_dDdT2 = 3.5988912573413091E-004
      real(pr) :: test_dDi(nc) = [45.043635428115536, 65.731805235540278]
      real(pr) :: test_dDidT(nc) = [-8.2934765998547377E-002, -0.18941780037882283]
      real(pr) :: test_dDij(nc, nc) = reshape( &
                  [25.390034959603035, 49.957034706240584, &
                   49.957034706240584, 69.675496643514933], [nc, nc])

      real(pr) :: ai(nc), daidt(nc), daidt2(nc)
      real(pr) :: n(nc), T, Tr(nc), Tc(nc)

      real(pr) :: D, dDdT, dDdT2, dDi(nc), dDidT(nc), dDij(nc,nc)
      
      type(CubicEoS) :: model

      model = binary_NRTL_SRK()

      n = [0.2, 0.8]
      T = 150
      Tc = model%components%Tc
      Tr = T/Tc

      call model%alpha%alpha(Tr, ai, daidt, daidt2)
      ai = ai*model%ac
      daidt = daidt*model%ac/Tc
      daidt2 = daidt2*model%ac/Tc**2

      call model%mixrule%Dmix(n, T, ai, daidt, daidt2, D, dDdT, dDdT2, dDi, dDidT, dDij)
      
      call check(error, allclose([D], [test_D], absolute_tolerance))
      call check(error, allclose([dDdT], [test_dDdT], absolute_tolerance))
      call check(error, allclose([dDdT2], [test_dDdT2], absolute_tolerance))
      call check(error, allclose([dDi], [test_dDi], absolute_tolerance))
      call check(error, allclose([dDidT], [test_dDidT], absolute_tolerance))
      call check(error, allclose([dDij], [test_dDij], absolute_tolerance))
   end subroutine

end module test_cubic_mixrules
