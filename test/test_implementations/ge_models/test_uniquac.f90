module test_uniquac
   use yaeos, only: pr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use auxiliar_functions, only: allclose, rel_error
   implicit none

contains
   subroutine collect_suite(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("Test UNIQUAC consistency mixture", test_uniquac_cons_mix), &
         new_unittest("Test UNIQUAC consistency pure", test_uniquac_cons_pure), &
         new_unittest("Test UNIQUAC against Caleb Bell's Thermo lib", test_against_caleb_thermo), &
         new_unittest("Test UNIQUAC against Gmehling et al. book", test_against_gmehling), &
         new_unittest("Test UNIQUAC against Gmehling et al. book 2", test_against_gmehling2), &
         new_unittest("Test UNIQUAC Temperature dependence", test_temperature_dependence) &
         ]
   end subroutine collect_suite

   subroutine test_uniquac_cons_mix(error)
      use yaeos, only: pr, R
      use yaeos, only: Groups, setup_uniquac, UNIQUAC
      use yaeos__consistency_gemodel, only: ge_consistency
      use yaeos__consistency_gemodel, only: numeric_ge_derivatives

      type(error_type), allocatable, intent(out) :: error

      type(UNIQUAC) :: model

      integer, parameter :: nc = 3

      real(pr) :: Ge, Gen(nc), GeT, GeT2, GeTn(nc), Gen2(nc, nc)
      real(pr) :: Gei, Geni(nc), GeTi, GeT2i, GeTni(nc), Gen2i(nc, nc)
      real(pr) :: Ge_n, Gen_n(nc), GeT_n, GeT2_n, GeTn_n(nc), Gen2_n(nc, nc)

      real(pr) :: n(nc), T, rs(nc), qs(nc)
      real(pr) :: A(nc,nc), B(nc,nc), C(nc,nc), D(nc,nc), E(nc,nc)
      real(pr) :: dt, dn

      real(pr) :: eq58, eq59(nc), eq60(nc,nc), eq61(nc)

      integer :: i, j

      T = 350.15_pr
      n = [5.0_pr, 8.0_pr, 10.0_pr]

      dt = 0.1_pr
      dn = 0.001_pr

      A(1,:) = [0.0_pr, -75.46_pr, -60.15_pr]
      A(2,:) = [120.20_pr, 0.0_pr, 44.22_pr]
      A(3,:) = [120.20_pr, 33.21_pr, 0.0_pr]

      B(1,:) = [0.0_pr, -0.10062_pr, 0.2566_pr]
      B(2,:) = [0.44835_pr, 0.0_pr, -0.01325_pr]
      B(3,:) = [0.44835_pr, 0.124_pr, 0.0_pr]

      C(1,:) = [0.0_pr, -0.0008052_pr, 0.00021_pr]
      C(2,:) = [0.0004704_pr, 0.0_pr, -0.00033_pr]
      C(3,:) = [0.0004704_pr, -0.000247_pr, 0.0_pr]

      D(1,:) = [0.0_pr, -0.001_pr, 0.0002_pr]
      D(2,:) = [-0.001_pr, 0.0_pr, 0.0002_pr]
      D(3,:) = [-0.001_pr, 0.0002_pr, 0.0_pr]

      E(1,:) = [0.0_pr, -0.00001_pr, 0.00001_pr]
      E(2,:) = [-0.00001_pr, 0.0_pr, 0.00001_pr]
      E(3,:) = [-0.00001_pr, 0.00001_pr, 0.0_pr]

      rs = [0.92_pr, 2.1055_pr, 1.5_pr]
      qs = [1.4_pr, 1.972_pr, 1.4_pr]

      model = setup_uniquac(qs, rs, A, B, C, D, E)

      ! ========================================================================
      ! Call analytic derivatives
      ! ------------------------------------------------------------------------
      call model%excess_gibbs(n, T, Ge=Ge, GeT=GeT, GeT2=GeT2, Gen=Gen, Gen2=Gen2, GeTn=GeTn)
      call model%excess_gibbs(n, T, Ge=Gei)
      call model%excess_gibbs(n, T, GeT=GeTi)
      call model%excess_gibbs(n, T, GeT2=GeT2i)
      call model%excess_gibbs(n, T, Gen=Geni)
      call model%excess_gibbs(n, T, GeTn=GeTni)
      call model%excess_gibbs(n, T, Gen2=Gen2i)

      ! ========================================================================
      ! Test single calls
      ! ------------------------------------------------------------------------
      call check(error, abs(Ge - Gei) < 1e-12_pr)
      call check(error, abs(GeT - GeTi) < 1e-12_pr)
      call check(error, abs(GeT2 - GeT2i) < 1e-12_pr)
      call check(error, allclose(Gen, Geni, 1e-12_pr))
      call check(error, allclose(GeTn, GeTni, 1e-12_pr))
      call check(error, allclose(Gen2(1,:), Gen2i(1,:), 1e-12_pr))
      call check(error, allclose(Gen2(2,:), Gen2i(2,:), 1e-12_pr))

      ! ========================================================================
      ! Test pair calls
      ! ------------------------------------------------------------------------
      ! Ge
      call model%excess_gibbs(n, T, Ge=Gei, GeT=GeTi)
      call check(error, abs(Ge - Gei) <= 1e-10)
      call check(error, abs(GeT - GeTi) <= 1e-10)

      call model%excess_gibbs(n, T, Ge=Gei, GeT2=GeT2i)
      call check(error, abs(Ge - Gei) <= 1e-10)
      call check(error, abs(GeT2 - GeT2i) <= 1e-10)

      call model%excess_gibbs(n, T, Ge=Gei, Gen=Geni)
      call check(error, abs(Ge - Gei) <= 1e-10)
      call check(error, allclose(Gen, Geni, 1e-10_pr))

      call model%excess_gibbs(n, T, Ge=Gei, GeTn=GeTni)
      call check(error, abs(Ge - Gei) <= 1e-10)
      call check(error, allclose(GeTn, GeTni, 1e-10_pr))

      call model%excess_gibbs(n, T, Ge=Gei, Gen2=Gen2i)
      call check(error, abs(Ge - Gei) <= 1e-10)
      call check(error, allclose(Gen2(1,:), Gen2i(1,:), 1e-10_pr))
      call check(error, allclose(Gen2(2,:), Gen2i(2,:), 1e-10_pr))

      ! Ge_T
      call model%excess_gibbs(n, T, GeT=GeTi, GeT2=GeT2i)
      call check(error, abs(GeT - GeTi) <= 1e-10)
      call check(error, abs(GeT2 - GeT2i) <= 1e-10)

      call model%excess_gibbs(n, T, GeT=GeTi, Gen=Geni)
      call check(error, abs(GeT - GeTi) <= 1e-10)
      call check(error, allclose(Gen, Geni, 1e-10_pr))

      call model%excess_gibbs(n, T, GeT=GeTi, GeTn=GeTni)
      call check(error, abs(GeT - GeTi) <= 1e-10)
      call check(error, allclose(GeTn, GeTni, 1e-10_pr))

      call model%excess_gibbs(n, T, GeT=GeTi, Gen2=Gen2i)
      call check(error, abs(GeT - GeTi) <= 1e-10)
      call check(error, allclose(Gen2(1,:), Gen2i(1,:), 1e-10_pr))
      call check(error, allclose(Gen2(2,:), Gen2i(2,:), 1e-10_pr))

      ! Ge_T2
      call model%excess_gibbs(n, T, GeT2=GeT2i, Gen=Geni)
      call check(error, abs(GeT2 - GeT2i) <= 1e-10)
      call check(error, allclose(Gen, Geni, 1e-10_pr))

      call model%excess_gibbs(n, T, GeT2=GeT2i, GeTn=GeTni)
      call check(error, abs(GeT2 - GeT2i) <= 1e-10)
      call check(error, allclose(GeTn, GeTni, 1e-10_pr))

      call model%excess_gibbs(n, T, GeT2=GeT2i, Gen2=Gen2i)
      call check(error, abs(GeT2 - GeT2i) <= 1e-10)
      call check(error, allclose(Gen2(1,:), Gen2i(1,:), 1e-10_pr))
      call check(error, allclose(Gen2(2,:), Gen2i(2,:), 1e-10_pr))

      ! Geni
      call model%excess_gibbs(n, T, Gen=Geni, GeTn=GeTni)
      call check(error, allclose(Gen, Geni, 1e-10_pr))
      call check(error, allclose(GeTn, GeTni, 1e-10_pr))

      call model%excess_gibbs(n, T, Gen=Geni, Gen2=Gen2i)
      call check(error, allclose(Gen, Geni, 1e-10_pr))
      call check(error, allclose(Gen2(1,:), Gen2i(1,:), 1e-10_pr))
      call check(error, allclose(Gen2(2,:), Gen2i(2,:), 1e-10_pr))

      ! ========================================================================
      ! Just one triplet call test
      ! ------------------------------------------------------------------------
      call model%excess_gibbs(n, T, Ge=Gei, GeT=GeTi, Gen2=Gen2i)
      call check(error, abs(Ge - Gei) <= 1e-10)
      call check(error, abs(GeT - GeTi) <= 1e-10)
      call check(error, allclose(Gen2(1,:), Gen2i(1,:), 1e-10_pr))
      call check(error, allclose(Gen2(2,:), Gen2i(2,:), 1e-10_pr))

      ! ========================================================================
      ! Call numeric derivatives
      ! ------------------------------------------------------------------------
      call numeric_ge_derivatives(model, n, T, dn, 0.01_pr, Ge=Ge_n, GeT=GeT_n)
      call numeric_ge_derivatives(model, n, T, dn, dt, Ge=Ge_n, Gen=Gen_n)
      call numeric_ge_derivatives(model, n, T, dn, dt, Ge=Ge_n, GeT2=GeT2_n)
      call numeric_ge_derivatives(model, n, T, dn, dt, Ge=Ge_n, GeTn=GeTn_n)
      call numeric_ge_derivatives(model, n, T, 0.01_pr, dt, Ge=Ge_n, Gen2=Gen2_n)

      ! Derivatives checks
      call check(error, abs(Ge - Ge_n) < 1e-10)
      call check(error, abs(GeT - GeT_n) < 1e-6)
      call check(error, allclose(Gen, Gen_n, 1e-6_pr))
      call check(error, abs(GeT2 - GeT2_n) < 1e-6)
      call check(error, allclose(GeTn, GeTn_n, 1e-6_pr))
      call check(error, allclose(Gen2(1,:), Gen2_n(1,:), 1e-5_pr))
      call check(error, allclose(Gen2(2,:), Gen2_n(2,:), 1e-5_pr))
      call check(error, allclose(Gen2(3,:), Gen2_n(3,:), 1e-5_pr))

      ! ========================================================================
      ! Consistency tests
      ! ------------------------------------------------------------------------
      call ge_consistency(model, n, t, eq58=eq58)
      call ge_consistency(model, n, t, eq59=eq59)
      call ge_consistency(model, n, t, eq60=eq60)
      call ge_consistency(model, n, t, eq61=eq61)

      ! Eq 58
      call check(error, abs(eq58) < 1e-10_pr)

      ! Eq 59
      do i=1,size(n)
         call check(error, abs(eq59(i)) < 1e-10_pr)
      end do

      ! Eq 60
      do i=1,size(n)
         do j=1,size(n)
            call check(error, abs(eq60(i, j)) < 1e-10_pr)
         end do
      end do

      ! Eq 61
      do i=1,size(n)
         call check(error, abs(eq61(i)) < 1e-10_pr)
      end do
   end subroutine test_uniquac_cons_mix

   subroutine test_uniquac_cons_pure(error)
      use yaeos, only: pr, UNIQUAC, setup_uniquac

      type(error_type), allocatable, intent(out) :: error

      type(UNIQUAC) :: model

      real(pr) :: Ge, Gen(1), GeT, GeT2, GeTn(1), Gen2(1, 1), ln_gammas(1)
      real(pr) :: T, n(1), rs(1), qs(1)

      T = 303.15
      n = [400.0]

      rs = [2.1055_pr]
      qs = [1.972_pr]

      model  = setup_uniquac(qs, rs)

      ! Evaluate Ge and derivatives
      call model%excess_gibbs(n, T, Ge, GeT, GeT2, Gen, GeTn, Gen2)
      call model%ln_activity_coefficient(n, T, ln_gammas)

      ! All must be zero for a pure compounds
      call check(error, abs(Ge) < 1e-10_pr)
      call check(error, abs(GeT) < 1e-10_pr)
      call check(error, abs(GeT2) < 1e-10_pr)

      call check(error, abs(Gen(1)) < 1e-10_pr)
      call check(error, abs(GeTn(1)) < 1e-10_pr)
      call check(error, abs(ln_gammas(1)) < 1e-10_pr)

      call check(error, abs(Gen2(1, 1)) < 1e-10_pr)
   end subroutine test_uniquac_cons_pure

   subroutine test_against_caleb_thermo(error)
      ! https://github.com/CalebBell/thermo
      use yaeos, only: pr, R
      use yaeos, only: Groups, setup_uniquac, UNIQUAC

      type(error_type), allocatable, intent(out) :: error

      type(UNIQUAC) :: model

      integer, parameter :: nc = 3

      real(pr) :: Ge, Gen(nc), GeT, GeT2, GeTn(nc), Gen2(nc, nc)
      real(pr) :: ln_gammas(nc)
      real(pr) :: rs(nc), qs(nc)
      real(pr), dimension(nc,nc) :: A, B, C, D, F, E

      real(pr) :: n(nc), T, n_t

      T = 298.15_pr
      n = [20.0_pr, 70.0_pr, 10.0_pr]
      n_t = sum(n)

      A(1,:) = [0.0_pr, -75.46_pr, -60.15_pr]
      A(2,:) = [120.20_pr, 0.0_pr, 44.22_pr]
      A(3,:) = [120.20_pr, 33.21_pr, 0.0_pr]

      B(1,:) = [0.0_pr, -0.10062_pr, 0.2566_pr]
      B(2,:) = [0.44835_pr, 0.0_pr, -0.01325_pr]
      B(3,:) = [0.44835_pr, 0.124_pr, 0.0_pr]

      C(1,:) = [0.0_pr, -0.0008052_pr, 0.00021_pr]
      C(2,:) = [0.0004704_pr, 0.0_pr, -0.00033_pr]
      C(3,:) = [0.0004704_pr, -0.000247_pr, 0.0_pr]

      D(1,:) = [0.0_pr, -0.001_pr, 0.0002_pr]
      D(2,:) = [-0.001_pr, 0.0_pr, 0.0002_pr]
      D(3,:) = [-0.001_pr, 0.0002_pr, 0.0_pr]

      E(1,:) = [0.0_pr, -0.00001_pr, 0.00001_pr]
      E(2,:) = [-0.00001_pr, 0.0_pr, 0.00001_pr]
      E(3,:) = [-0.00001_pr, 0.00001_pr, 0.0_pr]

      rs = [0.92_pr, 2.1055_pr, 1.5_pr]
      qs = [1.4_pr, 1.972_pr, 1.4_pr]

      ! setup UNIQUAC model
      model = setup_uniquac(qs, rs, A, B, C, D, E)

      ! Call all Ge and derivatives
      call model%excess_gibbs(n, T, Ge, GeT, GeT2, Gen, GeTn, Gen2)

      ! Call GeModel class method
      call model%ln_activity_coefficient(n, T, ln_gammas)

      ! ========================================================================
      ! Test against Caleb Bell's implementation
      ! ------------------------------------------------------------------------
      ! Ge
      call check(error, abs(Ge / n_t - (-2060.2989541519596_pr)) <= 1e-7)

      ! Gen
      call check(error, allclose(&
         Gen/R/T, [-164.62277497059728_pr, -60.906444787104235_pr, -75.52457152449654_pr], 1e-8_pr))

      ! ! ln_gammas
      call check(error, allclose(ln_gammas, [-164.62277497059728_pr, -60.906444787104235_pr, -75.52457152449654_pr], 1e-8_pr))

      ! Gen2
      call check(error, allclose(Gen2(1,:)*100, [1685.266356403747_pr, -484.7600670734814_pr, 22.787756706880828_pr]/n_t, 1e-7_pr))
      call check(error, allclose(Gen2(2,:)*100, [-484.76006707348006_pr, 7190.325877999089_pr, -49362.76101184667_pr]/n_t, 1e-7_pr))
      call check(error, allclose(Gen2(3,:)*100, [22.78775670689265_pr, -49362.76101184663_pr, 345493.75156951265_pr]/n_t, 1e-7_pr))

      ! GeT
      call check(error, abs(GeT * 100 / n_t - (-709.4126206919739_pr)) < 1e-7)

      ! GeT2
      call check(error, abs(GeT2*100/n_t - (-0.18488719888568747_pr)) < 1e-10)

      ! GeTn
      call check(error, allclose(GeTn*100, [-1344.5725109529376_pr, -536.5213350370494_pr, -649.331839754517_pr], 1e-8_pr))
   end subroutine test_against_caleb_thermo

   subroutine test_against_gmehling(error)
      ! n-butanol - water from Chemical Thermodynamics for Process Simulation
      ! Gmehling et al.
      use yaeos, only: pr, R
      use yaeos, only: setup_uniquac, UNIQUAC

      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 2

      real(pr) :: rs(nc), qs(nc)
      real(pr) :: b(nc, nc), b12, b21
      real(pr) :: x1(9), x2(9)
      real(pr) :: a1_g(9), a2_g(9)

      real(pr) :: ln_gammas(nc), T

      integer :: i

      type(UNIQUAC) :: model

      rs = [3.4543_pr, 0.92_pr]
      qs = [3.052_pr, 1.4_pr]

      T = 323.15_pr

      ! Calculate bij from DUij. We need -DU/R to get bij
      b12 = -129.7_pr * 0.04184_pr / R ! cal to bar L
      b21 = -489.6_pr * 0.04184_pr / R ! cal to bar L

      b(1,:) = [0.0_pr, b12]
      b(2,:) = [b21, 0.0_pr]

      model = setup_uniquac(qs, rs, bij=b)

      ! Test against Gmehling et al. book
      x1 = [0.005_pr, 0.01_pr, 0.015_pr, 0.02_pr, 0.05_pr, 0.1_pr, 0.2_pr, 0.4_pr, 0.6_pr]
      x2 = 1.0_pr - x1

      a1_g = [0.2824_pr, 0.4972_pr, 0.6598_pr, 0.9801_pr, 1.0495_pr, 0.9605_pr, 0.7253_pr, 0.6141_pr, 0.6802_pr]
      a2_g = [0.9953_pr, 0.9913_pr, 0.9878_pr, 0.9848_pr, 0.9766_pr, 0.9842_pr, 1.0332_pr, 1.0951_pr, 0.9766_pr]

      do i=1,9
         call model%ln_activity_coefficient([x1(i), x2(i)], T, ln_gammas)

         call check(error, abs(exp(ln_gammas(1))*x1(i) - a1_g(i)) / a1_g(i) < 0.2)
         call check(error, abs(exp(ln_gammas(2))*x2(i) - a2_g(i)) / a2_g(i) < 0.2)
      end do
   end subroutine test_against_gmehling

   subroutine test_against_gmehling2(error)
      ! water - ethanol- benzene from Chemical Thermodynamics for Process
      ! Simulation Gmehling et al.
      use yaeos, only: pr, R
      use yaeos, only: setup_uniquac, UNIQUAC

      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 3

      real(pr) :: rs(nc), qs(nc)
      real(pr) :: b(nc, nc)
      real(pr) :: z(nc)

      real(pr) :: ln_gammas(nc), T

      type(UNIQUAC) :: model

      rs = [0.92_pr, 2.1055_pr, 3.1878_pr]
      qs = [1.4_pr, 1.972_pr, 2.4_pr]

      T = 298.15_pr

      ! Calculate bij from DUij. We need -DU/R to get bij
      b(1,:) = [0.0_pr, -526.02_pr, -309.64_pr]
      b(2,:) = [318.06_pr, 0.0_pr, 91.532_pr]
      b(3,:) = [-1325.1_pr, -302.57_pr, 0.0_pr]

      model = setup_uniquac(qs, rs, bij=b)

      ! Test against Gmehling et al. book
      z = [0.8_pr, 0.1_pr, 0.2_pr]

      call model%ln_activity_coefficient(z, T, ln_gammas)

      call check(error, allclose(exp(ln_gammas), [1.570_pr, 0.2948_pr, 18.11_pr], 1e-3_pr))

      z = [0.2_pr, 0.2_pr, 0.8_pr]

      call model%ln_activity_coefficient(z, T, ln_gammas)

      call check(error, allclose(exp(ln_gammas), [8.856_pr, 0.860_pr, 1.425_pr], 1e-3_pr))

   end subroutine test_against_gmehling2

   subroutine test_temperature_dependence(error)
      use yaeos, only: pr, R
      use yaeos, only: setup_uniquac, UNIQUAC

      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 3

      real(pr) :: rs(nc), qs(nc)
      real(pr) :: n(nc)

      real(pr) :: ln_gammas_p(nc), Ge_p, GeT_p, GeT2_p, Gen_p(nc)
      real(pr) :: GeTn_p(nc), Gen2_p(nc, nc)

      real(pr) :: ln_gammas_t(nc), Ge_t, GeT_t, GeT2_t, Gen_t(nc)
      real(pr) :: GeTn_t(nc), Gen2_t(nc, nc)

      ! models
      type(UNIQUAC) :: model20
      type(UNIQUAC) :: model30
      type(UNIQUAC) :: model_const
      type(UNIQUAC) :: model_lineal
      type(UNIQUAC) :: model_ln
      type(UNIQUAC) :: model_quad

      ! Parameters
      real(pr) :: bij20(nc, nc), bij30(nc, nc), aij(nc, nc), dij(nc, nc)
      real(pr) :: cij(nc, nc), eij(nc, nc)

      n = [5.0_pr, 15.0_pr, 65.0_pr]
      rs = [1.972_pr, 10.496_pr, 31.764_pr]
      qs = [2.105_pr, 12.746_pr, 39.178_pr]

      ! =======================================================================
      ! Models with constants parameters to test
      ! -----------------------------------------------------------------------
      bij20(1,:) = [0.0_pr, 43.552026684128634_pr, 97.82405844415928_pr]
      bij20(2,:) = [-48.846029213395745_pr, 0.0_pr, 180.94738666155268_pr]
      bij20(3,:) = [-384.11635874542793_pr, -208.59497463051014_pr, 0.0_pr]

      model20 = setup_uniquac(qs, rs, bij=bij20)

      bij30(1,:) = [0.0_pr, 21.657910812761273_pr, 99.58458376639004_pr]
      bij30(2,:) = [-37.74765519959855_pr, 0.0_pr, 186.155606345583_pr]
      bij30(3,:) = [-379.32149369269956_pr, -233.34490510114676_pr, 0.0_pr]

      model30 = setup_uniquac(qs, rs, bij=bij30)

      ! =======================================================================
      ! Tests
      ! -----------------------------------------------------------------------
      ! Test against a non temperature dependent model
      ! 20 C
      model_const = setup_uniquac(qs, rs, aij=bij20 / 293.15_pr)

      call model_const%ln_activity_coefficient(n, 293.15_pr, ln_gammas_t)
      call model_const%excess_gibbs(&
         n, 293.15_pr, Ge_t, GeT_t, GeT2_t, Gen_t, GeTn_t, Gen2_t &
         )

      call model20%ln_activity_coefficient(n, 293.15_pr, ln_gammas_p)
      call model20%excess_gibbs(&
         n, 293.15_pr, Ge_p, GeT_p, GeT2_p, Gen_p, GeTn_p, Gen2_p &
         )

      call check(error, allclose(ln_gammas_t, ln_gammas_p, 1e-10_pr))
      call check(error, abs(Ge_t - Ge_p) < 1e-10_pr)
      call check(error, abs(GeT_t - Ge_t / 293.15_pr) < 1e-10_pr) ! funny
      call check(error, abs(GeT2_t) < 1e-10_pr) ! funny x2
      call check(error, allclose(Gen_t, Gen_p, 1e-10_pr))
      call check(error, allclose(GeTn_t, Gen_t / 293.15_pr, 1e-10_pr)) ! funny3
      call check(error, allclose(Gen2_t(1,:), Gen2_p(1,:), 1e-10_pr))
      call check(error, allclose(Gen2_t(2,:), Gen2_p(2,:), 1e-10_pr))
      call check(error, allclose(Gen2_t(3,:), Gen2_p(3,:), 1e-10_pr))

      ! 30 C
      model_const = setup_uniquac(qs, rs, aij=bij30 / 303.15_pr)

      call model_const%ln_activity_coefficient(n, 303.15_pr, ln_gammas_t)
      call model_const%excess_gibbs(&
         n, 303.15_pr, Ge_t, GeT_t, GeT2_t, Gen_t, GeTn_t, Gen2_t &
         )

      call model30%ln_activity_coefficient(n, 303.15_pr, ln_gammas_p)
      call model30%excess_gibbs(&
         n, 303.15_pr, Ge_p, GeT_p, GeT2_p, Gen_p, GeTn_p, Gen2_p &
         )

      call check(error, allclose(ln_gammas_t, ln_gammas_p, 1e-10_pr))
      call check(error, abs(Ge_t - Ge_p) < 1e-10_pr)
      call check(error, abs(GeT_t - Ge_t / 303.15_pr) < 1e-10_pr)
      call check(error, abs(GeT2_t) < 1e-10_pr)
      call check(error, allclose(Gen_t, Gen_p, 1e-10_pr))
      call check(error, allclose(GeTn_t, Gen_t / 303.15_pr, 1e-10_pr))
      call check(error, allclose(Gen2_t(1,:), Gen2_p(1,:), 1e-10_pr))
      call check(error, allclose(Gen2_t(2,:), Gen2_p(2,:), 1e-10_pr))
      call check(error, allclose(Gen2_t(3,:), Gen2_p(3,:), 1e-10_pr))

      ! Test against a linear temperature dependent model
      aij(1,:) = [0.0_pr, 2.4094201446651944_pr, 0.4861465075816882_pr]
      aij(2,:) = [-1.4009801734684584_pr, 0.0_pr, 0.710500847827416_pr]
      aij(3,:) = [-3.0410597123328746_pr, 0.9936949460465081_pr, 0.0_pr]

      dij(1,:) = [0.0_pr, -0.007712278604397001_pr, -0.0005200301448241018_pr]
      dij(2,:) = [0.004210661705262751_pr, 0.0_pr, -0.0003180930394505843_pr]
      dij(3,:) = [0.005903984936450968_pr, -0.005817018267835272_pr, 0.0_pr]

      model_lineal = setup_uniquac(qs, rs, aij=aij, dij=dij)

      ! 20 C
      call model_lineal%ln_activity_coefficient(n, 293.15_pr, ln_gammas_t)
      call model_lineal%excess_gibbs(&
         n, 293.15_pr, Ge_t, GeT_t, GeT2_t, Gen_t, GeTn_t, Gen2_t &
         )

      call model20%ln_activity_coefficient(n, 293.15_pr, ln_gammas_p)
      call model20%excess_gibbs(&
         n, 293.15_pr, Ge_p, GeT_p, GeT2_p, Gen_p, GeTn_p, Gen2_p &
         )

      call check(error, allclose(ln_gammas_t, ln_gammas_p, 1e-10_pr))
      call check(error, abs(Ge_t - Ge_p) < 1e-10_pr)
      call check(error, allclose(Gen_t, Gen_p, 1e-10_pr))
      call check(error, allclose(Gen2_t(1,:), Gen2_p(1,:), 1e-10_pr))
      call check(error, allclose(Gen2_t(2,:), Gen2_p(2,:), 1e-10_pr))
      call check(error, allclose(Gen2_t(3,:), Gen2_p(3,:), 1e-10_pr))

      ! 30 C
      call model_lineal%ln_activity_coefficient(n, 303.15_pr, ln_gammas_t)
      call model_lineal%excess_gibbs(&
         n, 303.15_pr, Ge_t, GeT_t, GeT2_t, Gen_t, GeTn_t, Gen2_t &
         )

      call model30%ln_activity_coefficient(n, 303.15_pr, ln_gammas_p)
      call model30%excess_gibbs(&
         n, 303.15_pr, Ge_p, GeT_p, GeT2_p, Gen_p, GeTn_p, Gen2_p &
         )

      call check(error, allclose(ln_gammas_t, ln_gammas_p, 1e-10_pr))
      call check(error, abs(Ge_t - Ge_p) < 1e-10_pr)
      call check(error, allclose(Gen_t, Gen_p, 1e-10_pr))
      call check(error, allclose(Gen2_t(1,:), Gen2_p(1,:), 1e-10_pr))
      call check(error, allclose(Gen2_t(2,:), Gen2_p(2,:), 1e-10_pr))
      call check(error, allclose(Gen2_t(3,:), Gen2_p(3,:), 1e-10_pr))

      ! Test against a logarithmic temperature dependent model
      aij(1,:) = [0.0_pr, 13.209596948268368_pr, 1.214390103981535_pr]
      aij(2,:) = [-7.297537236522662_pr, 0.0_pr, 1.155954291911136_pr]
      aij(3,:) = [-11.308925077456808_pr, 9.139773295263835_pr, 0.0_pr]

      cij(1,:) = [0.0_pr, -2.2992002904891815_pr, -0.1550324516753101_pr]
      cij(2,:) = [1.2552910900252325_pr, 0.0_pr, -0.09483054830130232_pr]
      cij(3,:) = [1.7601080792377783_pr, -1.73418139790263_pr, 0.0_pr]

      model_ln = setup_uniquac(qs, rs, aij=aij, cij=cij)

      ! 20 C
      call model_ln%ln_activity_coefficient(n, 293.15_pr, ln_gammas_t)
      call model_ln%excess_gibbs(&
         n, 293.15_pr, Ge_t, GeT_t, GeT2_t, Gen_t, GeTn_t, Gen2_t &
         )

      call model20%ln_activity_coefficient(n, 293.15_pr, ln_gammas_p)
      call model20%excess_gibbs(&
         n, 293.15_pr, Ge_p, GeT_p, GeT2_p, Gen_p, GeTn_p, Gen2_p &
         )

      call check(error, allclose(ln_gammas_t, ln_gammas_p, 1e-10_pr))
      call check(error, abs(Ge_t - Ge_p) < 1e-10_pr)
      call check(error, allclose(Gen_t, Gen_p, 1e-10_pr))
      call check(error, allclose(Gen2_t(1,:), Gen2_p(1,:), 1e-10_pr))
      call check(error, allclose(Gen2_t(2,:), Gen2_p(2,:), 1e-10_pr))
      call check(error, allclose(Gen2_t(3,:), Gen2_p(3,:), 1e-10_pr))

      ! 30 C
      call model_ln%ln_activity_coefficient(n, 303.15_pr, ln_gammas_t)
      call model_ln%excess_gibbs(&
         n, 303.15_pr, Ge_t, GeT_t, GeT2_t, Gen_t, GeTn_t, Gen2_t &
         )

      call model30%ln_activity_coefficient(n, 303.15_pr, ln_gammas_p)
      call model30%excess_gibbs(&
         n, 303.15_pr, Ge_p, GeT_p, GeT2_p, Gen_p, GeTn_p, Gen2_p &
         )

      call check(error, allclose(ln_gammas_t, ln_gammas_p, 1e-10_pr))
      call check(error, abs(Ge_t - Ge_p) < 1e-10_pr)
      call check(error, allclose(Gen_t, Gen_p, 1e-10_pr))
      call check(error, allclose(Gen2_t(1,:), Gen2_p(1,:), 1e-10_pr))
      call check(error, allclose(Gen2_t(2,:), Gen2_p(2,:), 1e-10_pr))
      call check(error, allclose(Gen2_t(3,:), Gen2_p(3,:), 1e-10_pr))

      ! Test against a quadratic temperature dependent model
      aij(1,:) = [0.0_pr, 1.2600355505795615_pr, 0.40864481611268855_pr]
      aij(2,:) = [-0.773452312613418_pr, 0.0_pr, 0.6630944640873461_pr]
      aij(3,:) = [-2.1611706837127835_pr, 0.12676682745485435_pr, 0.0_pr]

      eij(1,:) = [0.0_pr, -1.293355459399128e-05_pr, -8.720948261346668e-07_pr]
      eij(2,:) = [7.0613142801656066e-06_pr, 0.0_pr, -5.33444641037371e-07_pr]
      eij(3,:) = [9.901031253481415e-06_pr, -9.755187435578186e-06_pr, 0.0_pr]

      model_quad = setup_uniquac(qs, rs, aij=aij, eij=eij)

      ! 20 C
      call model_quad%ln_activity_coefficient(n, 293.15_pr, ln_gammas_t)
      call model_quad%excess_gibbs(&
         n, 293.15_pr, Ge_t, GeT_t, GeT2_t, Gen_t, GeTn_t, Gen2_t &
         )

      call model20%ln_activity_coefficient(n, 293.15_pr, ln_gammas_p)
      call model20%excess_gibbs(&
         n, 293.15_pr, Ge_p, GeT_p, GeT2_p, Gen_p, GeTn_p, Gen2_p &
         )

      call check(error, allclose(ln_gammas_t, ln_gammas_p, 1e-10_pr))
      call check(error, abs(Ge_t - Ge_p) < 1e-10_pr)
      call check(error, allclose(Gen_t, Gen_p, 1e-10_pr))
      call check(error, allclose(Gen2_t(1,:), Gen2_p(1,:), 1e-10_pr))
      call check(error, allclose(Gen2_t(2,:), Gen2_p(2,:), 1e-10_pr))
      call check(error, allclose(Gen2_t(3,:), Gen2_p(3,:), 1e-10_pr))

      ! 30 C
      call model_quad%ln_activity_coefficient(n, 303.15_pr, ln_gammas_t)
      call model_quad%excess_gibbs(&
         n, 303.15_pr, Ge_t, GeT_t, GeT2_t, Gen_t, GeTn_t, Gen2_t &
         )

      call model30%ln_activity_coefficient(n, 303.15_pr, ln_gammas_p)
      call model30%excess_gibbs(&
         n, 303.15_pr, Ge_p, GeT_p, GeT2_p, Gen_p, GeTn_p, Gen2_p &
         )

      call check(error, allclose(ln_gammas_t, ln_gammas_p, 1e-10_pr))
      call check(error, abs(Ge_t - Ge_p) < 1e-10_pr)
      call check(error, allclose(Gen_t, Gen_p, 1e-10_pr))
      call check(error, allclose(Gen2_t(1,:), Gen2_p(1,:), 1e-10_pr))
      call check(error, allclose(Gen2_t(2,:), Gen2_p(2,:), 1e-10_pr))
      call check(error, allclose(Gen2_t(3,:), Gen2_p(3,:), 1e-10_pr))
   end subroutine test_temperature_dependence
end module test_uniquac
