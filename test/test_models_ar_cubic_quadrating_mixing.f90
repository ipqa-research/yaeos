program test_cubic_mixrules
   use yaeos__constants, only: pr
   use auxiliar_functions, only: allclose
   use testing_aux, only: assert, test_title
   implicit none

   real(pr) :: absolute_tolerance = 1e-5_pr

   print *, test_title("CUBIC MIXING RULES")
   call test_QMR_RKPR
   call test_QMR_KTDEP
   call test_MHV_SRK
   print *, ""

contains

   subroutine test_QMR_RKPR
      use yaeos__constants, only: pr
      use yaeos__models_ar_cubic_quadratic_mixing, only: QMR

      type(QMR) :: mixrule

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

      call assert(allclose([D1], [D1_val], absolute_tolerance), "QMR D1")
      call assert(allclose([dD1i], [dD1i_val], absolute_tolerance), "QMR dD1i")
      call assert(allclose([dD1ij], [dD1ij_val], absolute_tolerance), "QMR dD1ij")
   end subroutine test_QMR_RKPR

   subroutine test_QMR_KTDEP
      use yaeos__models_ar_cubic_quadratic_mixing, only: QMRTD
      use yaeos, only: CubicEoS
      use fixtures_models, only: binary_PR78
      type(QMRTD) :: mixrule
      type(CubicEoS) :: model
      integer, parameter :: nc=2
      real(pr) :: k0(nc, nc), kinf(nc, nc), Tref(nc, nc), k02(nc, nc), k01(nc,nc)
      real(pr) :: a(nc), dadt(nc), dadt2(nc)
      real(pr) :: a_val(nc), dadt_val(nc), dadt2_val(nc)

      real(pr) :: T

      integer :: i

      a_val = [0.89240453619503446, 1.0876334517589699]
      dadt_val = [-0.3098030339083731, -0.488526759221918]
      dadt2_val = [0.17150004752436407, 0.41260112463979509]

      k0   = reshape([0, 1, 1, 0], [nc, nc])
      kinf = reshape([0, 2, 2, 0], [nc, nc])
      Tref = reshape([190, 310, 310, 310], [nc, nc])
      mixrule = QMRTD(k0=k0, k=kinf, Tref=Tref)

      model = binary_PR78()

      T = 250
      call model%alpha%alpha(T/model%components%Tc, a, dadt, dadt2)
      call mixrule%aij(T, a, dadt, dadt2, k0, kinf, tref)

      call assert(allclose(a, a_val, absolute_tolerance), "QMR KDTEP")
      call assert(allclose(dadt, dadt_val, absolute_tolerance), "QMR KDTEP")
      call assert(allclose(dadt2, dadt2_val, absolute_tolerance), "QMR KDTEP")
   end subroutine

   subroutine test_MHV_SRK
      use yaeos__constants, only: pr
      use yaeos__models, only: CubicEoS, NRTL, MHV
      use fixtures_models, only: binary_NRTL_SRK

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

      call assert(allclose([D], [test_D], absolute_tolerance), "MHV_SRK")
      call assert(allclose([dDdT], [test_dDdT], absolute_tolerance), "MHV_SRK")
      call assert(allclose([dDdT2], [test_dDdT2], absolute_tolerance), "MHV_SRK")
      call assert(allclose([dDi], [test_dDi], absolute_tolerance), "MHV_SRK")
      call assert(allclose([dDidT], [test_dDidT], absolute_tolerance), "MHV_SRK")
      call assert(allclose([dDij], [test_dDij], absolute_tolerance), "MHV_SRK")
   end subroutine

end program test_cubic_mixrules
