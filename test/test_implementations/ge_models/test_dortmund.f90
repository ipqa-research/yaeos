module test_dortmund
   use yaeos, only: pr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use auxiliar_functions, only: allclose, rel_error
   implicit none

contains
   subroutine collect_suite(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("Test dortmund consistency mixture", test_dortmund_cons_mix), &
         new_unittest("Test dortmund consistency pure", test_dortmund_cons_pure), &
         new_unittest("Test dortmund against Caleb Bell's Thermo lib", test_against_caleb_thermo) &
         ]
   end subroutine collect_suite

   subroutine test_dortmund_cons_mix(error)
      use yaeos, only: pr, R
      use yaeos, only: Groups, setup_dortmund, UNIFAC
      use yaeos__consistency_gemodel, only: ge_consistency
      use yaeos__consistency_gemodel, only: numeric_ge_derivatives
      use yaeos__models_ge_group_contribution_model_parameters, only: GeGCModelParameters
      use yaeos__models_ge_group_contribution_dortmund_parameters, only: DortmundParameters

      type(error_type), allocatable, intent(out) :: error

      type(UNIFAC) :: model

      integer, parameter :: nc = 4, ng = 4

      type(Groups) :: molecules(nc)

      type(GeGCModelParameters) :: parameters

      real(pr) :: Ge, Gen(nc), GeT, GeT2, GeTn(nc), Gen2(nc, nc)
      real(pr) :: Ge_n, Gen_n(nc), GeT_n, GeT2_n, GeTn_n(nc), Gen2_n(nc, nc)

      real(pr) :: n(nc), T
      real(pr) :: dt, dn

      real(pr) :: eq58, eq59(size(n)), eq60(size(n),size(n)), eq61(size(n))

      integer :: i, j

      T = 303.15
      n = [400.0, 100.0, 300.0, 200.0]

      dt = 0.1_pr
      dn = 0.01_pr

      molecules(1)%groups_ids = [16]
      molecules(1)%number_of_groups = [1]

      molecules(2)%groups_ids = [1, 2]
      molecules(2)%number_of_groups = [2, 4]

      molecules(3)%groups_ids = [1, 2, 14]
      molecules(3)%number_of_groups = [1, 1, 1]

      molecules(4)%groups_ids = [1, 7, 8, 78, 79]
      molecules(4)%number_of_groups = [2, 1, 1, 3, 1]

      parameters = DortmundParameters()

      model = setup_dortmund(molecules, parameters)

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
      call check(error, abs(GeT - GeT_n) < 1e-4)
      call check(error, allclose(Gen, Gen_n, 1e-4_pr))
      call check(error, abs(GeT2 - GeT2_n) < 1e-4_pr)
      call check(error, allclose(GeTn, GeTn_n, 1e-3_pr))
      call check(error, allclose(Gen2(1,:), Gen2_n(1,:), 1e-1_pr))
      call check(error, allclose(Gen2(2,:), Gen2_n(2,:), 1e-1_pr))
      call check(error, allclose(Gen2(3,:), Gen2_n(3,:), 1e-1_pr))


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
   end subroutine test_dortmund_cons_mix

   subroutine test_dortmund_cons_pure(error)
      use yaeos, only: pr
      use yaeos, only: Groups, UNIFAC, setup_dortmund

      type(error_type), allocatable, intent(out) :: error

      type(UNIFAC) :: model
      type(Groups) :: molecules(1)

      real(pr) :: Ge, Gen(1), GeT, GeT2, GeTn(1), Gen2(1, 1), ln_gammas(1)
      real(pr) :: T, n(1)

      integer :: i, j

      T = 303.15
      n = [400.0]

      ! Hexane [CH3, CH2]
      molecules(1)%groups_ids = [1, 2]
      molecules(1)%number_of_groups = [2, 4]

      model  = setup_dortmund(molecules)

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
   end subroutine test_dortmund_cons_pure

   subroutine test_against_caleb_thermo(error)
      ! https://github.com/CalebBell/thermo
      use yaeos, only: pr, R
      use yaeos, only: Groups, setup_dortmund, UNIFAC

      type(error_type), allocatable, intent(out) :: error

      type(UNIFAC) :: model

      integer, parameter :: nc = 3, ng = 4

      type(Groups) :: molecules(nc)

      real(pr) :: Ge, Gen(nc), GeT, GeT2, GeTn(nc), Gen2(nc, nc)
      real(pr) :: Ge_i, Gen_i(nc), GeT_i, GeT2_i, GeTn_i(nc), Gen2_i(nc, nc)
      real(pr) :: ln_gammas(nc)

      real(pr) :: n(nc), T, n_t

      T = 303.15_pr
      n = [2.0_pr, 5.0_pr, 3.0_pr]
      n_t = sum(n)

      molecules(1)%groups_ids = [1, 2]
      molecules(1)%number_of_groups = [2, 4]

      molecules(2)%groups_ids = [1, 2, 14]
      molecules(2)%number_of_groups = [1, 1, 1]

      molecules(3)%groups_ids = [1, 7, 8, 78, 79]
      molecules(3)%number_of_groups = [2, 1, 1, 3, 1]

      ! setup Dortmund model
      model = setup_dortmund(molecules)

      ! Call all Ge and derivatives
      call model%excess_gibbs(n, T, Ge, GeT, GeT2, Gen, GeTn, Gen2)

      ! Call Ge and derivatives individually
      call model%excess_gibbs(n, T, Ge_i)
      call model%excess_gibbs(n, T, GeT=GeT_i)
      call model%excess_gibbs(n, T, GeT2=GeT2_i)
      call model%excess_gibbs(n, T, Gen=Gen_i)
      call model%excess_gibbs(n, T, GeTn=GeTn_i)
      call model%excess_gibbs(n, T, Gen2=Gen2_i)

      ! Call GeModel class method
      call model%ln_activity_coefficient(n, T, ln_gammas)

      ! ========================================================================
      ! Test against Caleb Bell's implementation
      ! ------------------------------------------------------------------------
      ! Ge
      ! print *, Ge/n_t
      call check(error, abs(Ge / n_t - (14.432048256607143_pr)) <= 1e-5)

      ! ln_gammas
      ! print *, Gen/R/T, ln_gammas
      call check(error, allclose(Gen / R / T, [0.727831865699655_pr, 0.49653097358793297_pr, 0.595827311587003_pr], 1e-5_pr))
      call check(error, allclose(ln_gammas, [0.727831865699655_pr, 0.49653097358793297_pr, 0.595827311587003_pr], 1e-5_pr))

      ! Gen2
      ! print *, Gen2(1,:)
      ! print *, Gen2(2,:)
      ! print *, Gen2(3,:)
      call check(error, allclose(Gen2(1,:), [-27.396545333109035_pr, 22.045031798841627_pr, -18.47735610933041_pr] / n_t, 1e-5_pr))
      call check(error, allclose(Gen2(2,:), [22.04503179884173_pr, -22.7202581288098_pr, 23.170409015455117_pr] / n_t, 1e-5_pr))
      call check(error, allclose(Gen2(3,:), [-18.47735610933015_pr, 23.17040901545502_pr, -26.299110952871235_pr] / n_t, 1e-5_pr))

      ! GeT
      ! print *, GeT/n_t
      call check(error, abs(GeT / n_t - 0.019210544360072263_pr) < 1e-5)

      ! GeT2
      ! print *, GeT2 / n_t
      call check(error, abs(GeT2 / n_t - (-0.000614035558615488_pr)) < 1e-5)

      ! GeTn
      ! print *, GeTn
      call check(error, allclose(GeTn, [0.023995407535622692_pr, 0.028259971730944012_pr, 0.0009382566249192064_pr], 1e-5_pr))

      ! ========================================================================
      ! Test individual calls
      ! ------------------------------------------------------------------------
      call check(error, abs(Ge - Ge_i) <= 1e-10)
      call check(error, abs(GeT - GeT_i) <= 1e-10)
      call check(error, abs(GeT2 - GeT2_i) <= 1e-10)
      call check(error, allclose(Gen, Gen_i, 1e-10_pr))
      call check(error, allclose(GeTn, GeTn_i, 1e-10_pr))

      call check(error, allclose(Gen2(1,:), Gen2_i(1,:), 1e-10_pr))
      call check(error, allclose(Gen2(2,:), Gen2_i(2,:), 1e-10_pr))
      call check(error, allclose(Gen2(3,:), Gen2_i(3,:), 1e-10_pr))

      ! ========================================================================
      ! Test pair calls
      ! ------------------------------------------------------------------------
      ! Ge
      call model%excess_gibbs(n, T, Ge=Ge_i, GeT=GeT_i)
      call check(error, abs(Ge - Ge_i) <= 1e-10)
      call check(error, abs(GeT - GeT_i) <= 1e-10)

      call model%excess_gibbs(n, T, Ge=Ge_i, GeT2=GeT2_i)
      call check(error, abs(Ge - Ge_i) <= 1e-10)
      call check(error, abs(GeT2 - GeT2_i) <= 1e-10)

      call model%excess_gibbs(n, T, Ge=Ge_i, Gen=Gen_i)
      call check(error, abs(Ge - Ge_i) <= 1e-10)
      call check(error, allclose(Gen, Gen_i, 1e-10_pr))

      call model%excess_gibbs(n, T, Ge=Ge_i, GeTn=GeTn_i)
      call check(error, abs(Ge - Ge_i) <= 1e-10)
      call check(error, allclose(GeTn, GeTn_i, 1e-10_pr))

      call model%excess_gibbs(n, T, Ge=Ge_i, Gen2=Gen2_i)
      call check(error, abs(Ge - Ge_i) <= 1e-10)
      call check(error, allclose(Gen2(1,:), Gen2_i(1,:), 1e-10_pr))
      call check(error, allclose(Gen2(2,:), Gen2_i(2,:), 1e-10_pr))
      call check(error, allclose(Gen2(3,:), Gen2_i(3,:), 1e-10_pr))

      ! Ge_T
      call model%excess_gibbs(n, T, GeT=GeT_i, GeT2=GeT2_i)
      call check(error, abs(GeT - GeT_i) <= 1e-10)
      call check(error, abs(GeT2 - GeT2_i) <= 1e-10)

      call model%excess_gibbs(n, T, GeT=GeT_i, Gen=Gen_i)
      call check(error, abs(GeT - GeT_i) <= 1e-10)
      call check(error, allclose(Gen, Gen_i, 1e-10_pr))

      call model%excess_gibbs(n, T, GeT=GeT_i, GeTn=GeTn_i)
      call check(error, abs(GeT - GeT_i) <= 1e-10)
      call check(error, allclose(GeTn, GeTn_i, 1e-10_pr))

      call model%excess_gibbs(n, T, GeT=GeT_i, Gen2=Gen2_i)
      call check(error, abs(GeT - GeT_i) <= 1e-10)
      call check(error, allclose(Gen2(1,:), Gen2_i(1,:), 1e-10_pr))
      call check(error, allclose(Gen2(2,:), Gen2_i(2,:), 1e-10_pr))
      call check(error, allclose(Gen2(3,:), Gen2_i(3,:), 1e-10_pr))

      ! Ge_T2
      call model%excess_gibbs(n, T, GeT2=GeT2_i, Gen=Gen_i)
      call check(error, abs(GeT2 - GeT2_i) <= 1e-10)
      call check(error, allclose(Gen, Gen_i, 1e-10_pr))

      call model%excess_gibbs(n, T, GeT2=GeT2_i, GeTn=GeTn_i)
      call check(error, abs(GeT2 - GeT2_i) <= 1e-10)
      call check(error, allclose(GeTn, GeTn_i, 1e-10_pr))

      call model%excess_gibbs(n, T, GeT2=GeT2_i, Gen2=Gen2_i)
      call check(error, abs(GeT2 - GeT2_i) <= 1e-10)
      call check(error, allclose(Gen2(1,:), Gen2_i(1,:), 1e-10_pr))
      call check(error, allclose(Gen2(2,:), Gen2_i(2,:), 1e-10_pr))
      call check(error, allclose(Gen2(3,:), Gen2_i(3,:), 1e-10_pr))

      ! Gen_i
      call model%excess_gibbs(n, T, Gen=Gen_i, GeTn=GeTn_i)
      call check(error, allclose(Gen, Gen_i, 1e-10_pr))
      call check(error, allclose(GeTn, GeTn_i, 1e-10_pr))

      call model%excess_gibbs(n, T, Gen=Gen_i, Gen2=Gen2_i)
      call check(error, allclose(Gen, Gen_i, 1e-10_pr))
      call check(error, allclose(Gen2(1,:), Gen2_i(1,:), 1e-10_pr))
      call check(error, allclose(Gen2(2,:), Gen2_i(2,:), 1e-10_pr))
      call check(error, allclose(Gen2(3,:), Gen2_i(3,:), 1e-10_pr))

      ! ========================================================================
      ! Just one triplet call test
      ! ------------------------------------------------------------------------
      call model%excess_gibbs(n, T, Ge=Ge_i, GeT=GeT_i, Gen2=Gen2_i)
      call check(error, abs(Ge - Ge_i) <= 1e-10)
      call check(error, abs(GeT - GeT_i) <= 1e-10)
      call check(error, allclose(Gen2(1,:), Gen2_i(1,:), 1e-10_pr))
      call check(error, allclose(Gen2(2,:), Gen2_i(2,:), 1e-10_pr))
      call check(error, allclose(Gen2(3,:), Gen2_i(3,:), 1e-10_pr))
   end subroutine test_against_caleb_thermo
end module test_dortmund
