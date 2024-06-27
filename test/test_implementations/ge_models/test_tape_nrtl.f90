module test_tape_nrtl
   use yaeos__constants, only: pr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use auxiliar_functions, only: allclose
   implicit none

   real(pr) :: absolute_tolerance = 1e-4_pr

contains

   subroutine collect_suite(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("Test NRTL consistency mixture", test_cons_nrtl_mix), &
         new_unittest("Test NRTL consistency pure", test_cons_nrtl_pure), &
         new_unittest("NRTL", test_nrtl) &
         ]
   end subroutine collect_suite

   subroutine test_cons_nrtl_mix(error)
      use yaeos, only: pr, NRTL
      use yaeos__consistency_gemodel, only: ge_consistency
      use yaeos__consistency_gemodel, only: numeric_ge_derivatives

      type(error_type), allocatable, intent(out) :: error

      type(NRTL) :: model

      real(pr) :: Ge, Gen(2), GeT, GeT2, GeTn(2), Gen2(2, 2)
      real(pr) :: Ge_n, Gen_n(2), GeT_n, GeT2_n, GeTn_n(2), Gen2_n(2, 2)
      real(pr) :: ln_gammas(2)

      real(pr) :: eq58, eq59(2), eq60(2,2), eq61(2)

      real(pr) :: n(2), T
      real(pr) :: dt, dn

      real(pr) :: a(2, 2), b(2, 2), c(2, 2)

      integer :: i, j

      T = 303.15
      n = [400.0, 100.0]

      dt = 0.1_pr
      dn = 0.1_pr

      a = 0; b = 0; c = 0

      a(1, 2) = 3.458
      a(2, 1) = -0.801

      b(1, 2) = -586.1
      b(2, 1) = 246.2

      c(1, 2) = 0.3
      c(2, 1) = 0.3

      model = NRTL(a, b, c)

      ! ========================================================================
      ! Call analytic derivatives
      ! ------------------------------------------------------------------------
      call model%excess_gibbs(n, T, Ge, GeT, GeT2, Gen, GeTn, Gen2)

      ! ========================================================================
      ! Call numeric derivatives
      ! ------------------------------------------------------------------------
      call numeric_ge_derivatives(model, n, T, dn, 0.01_pr, Ge=Ge_n, GeT=GeT_n)
      call numeric_ge_derivatives(model, n, T, dn, dt, Ge=Ge_n, Gen=Gen_n)
      call numeric_ge_derivatives(model, n, T, dn, dt, Ge=Ge_n, GeT2=GeT2_n)
      call numeric_ge_derivatives(model, n, T, dn, dt, Ge=Ge_n, GeTn=GeTn_n)
      call numeric_ge_derivatives(model, n, T, dn, dt, Ge=Ge_n, Gen2=Gen2_n)

      ! Derivatives checks
      call check(error, abs(Ge - Ge_n) < 1e-10)
      call check(error, abs(GeT - GeT_n) < 1e-6)
      call check(error, allclose(Gen, Gen_n, 1e-6_pr))
      call check(error, abs(GeT2 - GeT2_n) < 1e-6)
      call check(error, allclose(GeTn, GeTn_n, 1e-6_pr))
      call check(error, allclose(Gen2(1,:), Gen2_n(1,:), 1e-5_pr))
      call check(error, allclose(Gen2(2,:), Gen2_n(2,:), 1e-5_pr))

      ! ========================================================================
      ! Consistency tests
      ! ------------------------------------------------------------------------
      call ge_consistency(&
         model, n, t, eq58=eq58, eq59=eq59, eq60=eq60, eq61=eq61 &
         )

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
   end subroutine test_cons_nrtl_mix

   subroutine test_cons_nrtl_pure(error)
      use yaeos, only: pr, NRTL

      type(error_type), allocatable, intent(out) :: error

      type(NRTL) :: model

      real(pr) :: Ge, Gen(1), GeT, GeT2, GeTn(1), Gen2(1, 1), ln_gammas(1)

      real(pr) :: n(1), T
      real(pr) :: dt, dn

      real(pr) :: a(1, 1), b(1, 1), c(1, 1)

      integer :: i, j

      T = 303.15
      n = [400.0]

      dt = 0.1_pr
      dn = 0.1_pr

      a = 0; b = 0; c = 0

      model = NRTL(a, b, c)

      ! Evaluate Ge and derivatives
      call model%excess_gibbs(n, T, Ge, GeT, GeT2, Gen, GeTn, Gen2)
      call model%ln_activity_coefficient(n, T, ln_gammas)

      ! All must be zero for a pure compounds
      call check(error, abs(Ge) < 1e-10_pr)
      call check(error, abs(GeT) < 1e-10_pr)
      call check(error, abs(GeT2) < 1e-10_pr)

      do i=1,size(n)
         call check(error, abs(Gen(i)) < 1e-10_pr)
         call check(error, abs(GeTn(i)) < 1e-10_pr)
         call check(error, abs(ln_gammas(i)) < 1e-10_pr)
      end do

      do i=1,size(n)
         do j=1,size(n)
            call check(error, abs(Gen2(i, j)) < 1e-10_pr)
         end do
      end do
   end subroutine test_cons_nrtl_pure

   subroutine test_nrtl(error)
      use yaeos__constants, only: pr
      use yaeos, only: GeModel, NRTL

      type(error_type), allocatable, intent(out) :: error

      type(NRTL) :: model

      integer, parameter :: n_c = 2

      real(pr) :: a(n_c, n_c), b(n_c, n_c), c(n_c, n_c)

      real(pr) :: n(n_c), T, n_t
      real(pr) :: Ge, GeT, GeT2
      real(pr) :: Gen(n_c), GeTn(n_c), Gen2(n_c, n_c)

      real(pr) :: Ge_val, GeT_val, GeT2_val
      real(pr) :: Gen_val(n_c), GeTn_val(n_c), Gen2_val(n_c**2)

      Ge_val = 0.73578738104034147
      GeT_val = 6.2478507144254576E-002
      GeT2_val = -8.6660745337171748E-004
      Gen_val = [1.4672850471801369, 0.42228836346995569]
      GeTn_val = [0.13831602133231552, 2.9976713504363872E-002]
      Gen2_val = [&
         -3.9849140290093517, 1.7078203950935533, &
         1.7078203950935529, -0.73192306801724871 &
         ]

      n = [3.0, 7.0]
      n_t = sum(n)
      T = 150

      a = 0; b = 0; c = 0

      a(1, 2) = 3.458
      a(2, 1) = -0.801

      b(1, 2) = -586.1
      b(2, 1) = 246.2

      c(1, 2) = 0.3
      c(2, 1) = 0.3

      model = NRTL(a, b, c)

      call model%excess_gibbs( &
         n, T, Ge=Ge, GeT=GeT, &
         GeT2=GeT2, Gen=Gen, GeTn=GeTn, Gen2=Gen2 &
         )

      call check(error, allclose([Ge / n_t], [Ge_val], absolute_tolerance))
      call check(error, allclose([GeT / n_t], [GeT_val], absolute_tolerance))
      call check(error, allclose([GeT2 / n_t], [GeT2_val], absolute_tolerance))
      call check(error, allclose([GeTn], [GeTn_val], absolute_tolerance))
      call check(error, allclose(Gen, Gen_val, absolute_tolerance))
      call check(error, allclose([Gen2], [Gen2_val / n_t], absolute_tolerance))
   end subroutine test_nrtl
end module test_tape_nrtl

