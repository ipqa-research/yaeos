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
         new_unittest("Test UNIQUAC against Gmehling et al. book 2", test_against_gmehling2) &
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
end module test_uniquac
