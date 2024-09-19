module test_psrk
   use yaeos, only: pr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use auxiliar_functions, only: allclose, rel_error
   implicit none

contains
   subroutine collect_suite(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("Test PSRK consistency mixture", test_psrk_cons_mix), &
         new_unittest("Test PSRK consistency pure", test_psrk_cons_pure), &
         new_unittest("Test PSRK against Caleb Bell's Thermo lib", test_against_caleb_thermo) &
         ]
   end subroutine collect_suite

   subroutine test_psrk_cons_mix(error)
      use yaeos, only: pr, R
      use yaeos, only: Groups, setup_psrk, UNIFAC
      use yaeos__consistency_gemodel, only: ge_consistency
      use yaeos__consistency_gemodel, only: numeric_ge_derivatives
      use yaeos__models_ge_group_contribution_model_parameters, only: GeGCModelParameters
      use yaeos__models_ge_group_contribution_psrk_parameters, only: PSRKParameters

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
      dn = 0.001_pr

      ! ! Hexane [CH3, CH2]
      molecules(1)%groups_ids = [1, 2]
      molecules(1)%number_of_groups = [2, 4]

      ! ! Ethanol [CH3, CH2, OH]
      molecules(2)%groups_ids = [1, 2, 14]
      molecules(2)%number_of_groups = [1, 1, 1]

      ! ! Toluene [ACH, ACCH3]
      molecules(3)%groups_ids = [9, 11]
      molecules(3)%number_of_groups = [5, 1]

      ! ! Cyclohexane [CH2]
      molecules(4)%groups_ids = [2]
      molecules(4)%number_of_groups = [6]

      parameters = PSRKParameters()

      model = setup_psrk(molecules, parameters)

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
      ! call check(error, abs(Ge - Ge_n) < 1e-10)
      ! print *, "Ge"
      ! print *, Ge, Ge_n
      ! call check(error, abs(GeT - GeT_n) < 1e-6)
      ! print *, "GeT"
      ! print *, GeT, Get_n
      ! call check(error, allclose(Gen, Gen_n, 1e-6_pr))
      ! print *, "Gen"
      ! print *, Gen
      ! print *, Gen_n
      ! call check(error, abs(GeT2 - GeT2_n) < 1e-6)
      ! print *, "GeT2"
      ! print *, GeT2, GeT2_n
      ! call check(error, allclose(GeTn, GeTn_n, 1e-6_pr))
      ! print *, "GeTn"
      ! print *, Getn
      ! print *, Getn_n
      
      ! print *, "Gen2"
      ! do i=1,3
      !    print *, GEn2(i, :)
      !    print *, GEn2_n(i, :)
      ! end do
      ! call check(error, allclose(Gen2(1,:), Gen2_n(1,:), 1e-4_pr))
      ! call check(error, allclose(Gen2(2,:), Gen2_n(2,:), 1e-4_pr))
      ! call check(error, allclose(Gen2(3,:), Gen2_n(3,:), 1e-4_pr))


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
   end subroutine test_psrk_cons_mix

   subroutine test_psrk_cons_pure(error)
      use yaeos, only: pr
      use yaeos, only: Groups, UNIFAC, setup_psrk

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

      model  = setup_psrk(molecules)

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
   end subroutine test_psrk_cons_pure

   subroutine test_against_caleb_thermo(error)
      ! https://github.com/CalebBell/thermo
      use yaeos, only: pr, R
      use yaeos, only: Groups, setup_psrk, UNIFAC

      type(error_type), allocatable, intent(out) :: error

      type(UNIFAC) :: model

      integer, parameter :: nc = 3, ng = 4

      type(Groups) :: molecules(nc)

      real(pr) :: Ge, Gen(nc), GeT, GeT2, GeTn(nc), Gen2(nc, nc)
      real(pr) :: Ge_i, Gen_i(nc), GeT_i, GeT2_i, GeTn_i(nc), Gen2_i(nc, nc)
      real(pr) :: ln_gammas(nc), psis(ng, ng)
      integer :: i

      real(pr) :: n(nc), T, n_t

      T = 150
      n = [20.0, 70.0, 10.0]
      n_t = sum(n)

      ! ! Ethane [CH3]
      molecules(1)%groups_ids = [1]
      molecules(1)%number_of_groups = [2]

      ! ! Ethanol [CH3, CH2, OH]
      molecules(2)%groups_ids = [1, 2, 14]
      molecules(2)%number_of_groups = [1, 1, 1]

      ! ! Methylamine [H3C-NH2]
      molecules(3)%groups_ids = [28]
      molecules(3)%number_of_groups = [1]

      ! setup PSRK model
      model = setup_psrk(molecules)

      ! Call all Ge and derivatives
      print *, "AAAAAAAAAA"
      call model%excess_gibbs(n, T, Ge, GeT, GeT2, Gen, GeTn, Gen2)
      call model%psi_function%psi(model%groups_stew, T, psi=psis)

      do i=1,ng
         print *, psis(i, :)
      end do
      call exit
      print *, "AAAAAAAAAA"

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

      print *, Ge/n_t
      call check(error, abs(Ge / n_t - (-3.223992676822129_pr)) <= 1e-10)

      ! Gen
      print *, Gen
      call check(error, allclose(Gen, [10.53032277_pr,  -2.37758326_pr, -36.65748951_pr], 1e-8_pr))

      ! ln_gammas
      call check(error, allclose(Gen / R / T, [0.84433781_pr, -0.19063836_pr, -2.93925506_pr], 1e-8_pr))
      call check(error, allclose(ln_gammas, [0.84433781_pr, -0.19063836_pr, -2.93925506_pr], 1e-8_pr))

      ! Gen2
      call check(error, allclose(Gen2(1,:), [-0.75249927_pr,  0.13440904_pr,  0.56413529_pr] * R*T/n_t, 1e-7_pr))
      call check(error, allclose(Gen2(2,:), [ 0.13440904_pr,  0.34708386_pr, -2.69840507_pr] * R*T/n_t, 1e-7_pr))
      call check(error, allclose(Gen2(3,:), [ 0.56413529_pr, -2.69840507_pr, 17.76056492_pr] * R*T/n_t, 1e-7_pr))

      ! dln_gammas_dn
      call check(error, allclose(Gen2(1,:) /R/T, [-0.752499273_pr, 0.134409037_pr, 0.564135287_pr] / n_t, 1e-7_pr))
      call check(error, allclose(Gen2(2,:) /R/T, [0.1344090359_pr, 0.3470838559_pr, -2.698405064_pr] / n_t, 1e-7_pr))
      call check(error, allclose(Gen2(3,:) /R/T, [0.5641352889_pr, -2.698405071_pr, 17.760564919_pr] / n_t, 1e-7_pr))

      ! GeT
      call check(error, abs(GeT / n_t - 0.03268447167877294_pr) < 1e-10)

      ! GeT2
      call check(error, abs(GeT2 / n_t - (-0.0003594405355829625_pr)) < 1e-10)

      ! GeTn
      call check(error, allclose(GeTn, [0.06015389_pr, 0.02239722_pr, 0.04975642_pr], 1e-6_pr))

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
end module test_psrk
