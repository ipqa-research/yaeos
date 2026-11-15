program test_nrtl
   use yaeos__constants, only: pr
   use auxiliar_functions, only: allclose
   use testing_aux, only: assert, test_title
   implicit none

   real(pr) :: absolute_tolerance = 1e-4_pr

   write(*, *) test_title("NRTL")
   call test_nrtl_individual_ge_calls()
   call test_cons_nrtl_mix()
   call test_cons_nrtl_pure()
   call test_nrtl_values()

contains

   subroutine test_nrtl_individual_ge_calls()
      use yaeos, only: pr, R
      use yaeos, only: NRTL
      use yaeos__consistency_gemodel, only: individual_ge_calls

      type(NRTL) :: model


      integer, parameter :: nc=2
      logical :: passed

      real(pr) :: n(nc), T, rs(nc), qs(nc)
      real(pr) :: a(nc,nc), b(nc,nc), c(nc,nc)

      T = 303.15
      n = [400.0, 100.0]

      a = 0; b = 0; c = 0

      a(1, 2) = 3.458
      a(2, 1) = -0.801

      b(1, 2) = -586.1
      b(2, 1) = 246.2

      c(1, 2) = 0.3
      c(2, 1) = 0.3

      model = NRTL(a, b, c)

      call individual_ge_calls(model, n, T, passed)

      call assert(passed, &
         "NRTL individual ge calls: passed")
   end subroutine test_nrtl_individual_ge_calls


   subroutine test_cons_nrtl_mix()
      use yaeos, only: pr, NRTL
      use yaeos__consistency_gemodel, only: ge_consistency
      use yaeos__consistency_gemodel, only: numeric_ge_derivatives


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
      call assert(abs(Ge - Ge_n) < 1e-10, &
         "NRTL consistency mixture: abs(Ge - Ge_n) < 1e-10")
      call assert(abs(GeT - GeT_n) < 1e-6, &
         "NRTL consistency mixture: abs(GeT - GeT_n) < 1e-6")
      call assert(allclose(Gen, Gen_n, 1e-6_pr), &
         "NRTL consistency mixture: allclose(Gen, Gen_n, 1e-6_pr)")
      call assert(abs(GeT2 - GeT2_n) < 1e-6, &
         "NRTL consistency mixture: abs(GeT2 - GeT2_n) < 1e-6")
      call assert(allclose(GeTn, GeTn_n, 1e-6_pr), &
         "NRTL consistency mixture: allclose(GeTn, GeTn_n, 1e-6_pr)")
      call assert(allclose(Gen2(1,:), Gen2_n(1,:), 1e-5_pr), &
         "NRTL consistency mixture: allclose(Gen2(1,:), Gen2_n(1,:), 1e-5_pr)")
      call assert(allclose(Gen2(2,:), Gen2_n(2,:), 1e-5_pr), &
         "NRTL consistency mixture: allclose(Gen2(2,:), Gen2_n(2,:), 1e-5_pr)")

      ! ========================================================================
      ! Consistency tests
      ! ------------------------------------------------------------------------
      call ge_consistency(&
         model, n, t, eq58=eq58, eq59=eq59, eq60=eq60, eq61=eq61 &
         )

      ! Eq 58
      call assert(abs(eq58) < 1e-10_pr, &
         "NRTL consistency mixture: abs(eq58) < 1e-10_pr")

      ! Eq 59
      do i=1,size(n)
         call assert(abs(eq59(i)) < 1e-10_pr, &
            "NRTL consistency mixture: abs(eq59(i)) < 1e-10_pr")
      end do

      ! Eq 60
      do i=1,size(n)
         do j=1,size(n)
            call assert(abs(eq60(i, j)) < 1e-10_pr, &
               "NRTL consistency mixture: abs(eq60(i, j)) < 1e-10_pr")
         end do
      end do

      ! Eq 61
      do i=1,size(n)
         call assert(abs(eq61(i)) < 1e-10_pr, &
            "NRTL consistency mixture: abs(eq61(i)) < 1e-10_pr")
      end do
   end subroutine test_cons_nrtl_mix

   subroutine test_cons_nrtl_pure()
      use yaeos, only: pr, NRTL


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
      call assert(abs(Ge) < 1e-10_pr, &
         "NRTL consistency pure: abs(Ge) < 1e-10_pr")
      call assert(abs(GeT) < 1e-10_pr, &
         "NRTL consistency pure: abs(GeT) < 1e-10_pr")
      call assert(abs(GeT2) < 1e-10_pr, &
         "NRTL consistency pure: abs(GeT2) < 1e-10_pr")

      do i=1,size(n)
         call assert(abs(Gen(i)) < 1e-10_pr, &
            "NRTL consistency pure: abs(Gen(i)) < 1e-10_pr")
         call assert(abs(GeTn(i)) < 1e-10_pr, &
            "NRTL consistency pure: abs(GeTn(i)) < 1e-10_pr")
         call assert(abs(ln_gammas(i)) < 1e-10_pr, &
            "NRTL consistency pure: abs(ln_gammas(i)) < 1e-10_pr")
      end do

      do i=1,size(n)
         do j=1,size(n)
            call assert(abs(Gen2(i, j)) < 1e-10_pr, &
               "NRTL consistency pure: abs(Gen2(i, j)) < 1e-10_pr")
         end do
      end do
   end subroutine test_cons_nrtl_pure

   subroutine test_nrtl_values()
      use yaeos__constants, only: pr
      use yaeos, only: GeModel, NRTL


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

      call assert(allclose([Ge / n_t], [Ge_val], absolute_tolerance), &
         "NRTL: allclose(ref, ref, absolute_tolerance)")
      call assert(allclose([GeT / n_t], [GeT_val], absolute_tolerance), &
         "NRTL: allclose(ref, ref, absolute_tolerance)")
      call assert(allclose([GeT2 / n_t], [GeT2_val], absolute_tolerance), &
         "NRTL: allclose(ref, ref, absolute_tolerance)")
      call assert(allclose([GeTn], [GeTn_val], absolute_tolerance), &
         "NRTL: allclose(ref, ref, absolute_tolerance)")
      call assert(allclose(Gen, Gen_val, absolute_tolerance), &
         "NRTL: allclose(Gen, Gen_val, absolute_tolerance)")
      call assert(allclose([Gen2], [Gen2_val / n_t], absolute_tolerance), &
         "NRTL: allclose(ref, ref, absolute_tolerance)")
   end subroutine test_nrtl_values
end program test_nrtl

