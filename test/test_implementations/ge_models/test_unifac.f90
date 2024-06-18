module test_unifac
   use yaeos, only: pr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use auxiliar_functions, only: allclose, rel_error
   implicit none

contains
   subroutine collect_suite(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("Test UNIFAC against Caleb Bell's Thermo lib", test_against_caleb_thermo) &
         ]
   end subroutine collect_suite

   subroutine test_against_caleb_thermo(error)
      ! https://github.com/CalebBell/thermo
      use yaeos, only: pr, R
      use yaeos, only: excess_gibbs, Groups, setup_unifac, UNIFAC

      type(error_type), allocatable, intent(out) :: error

      type(UNIFAC) :: model

      integer, parameter :: nc = 3, ng = 4

      type(Groups) :: molecules(nc)

      real(pr) :: Ge, Gen(nc), GeT, GeT2, GeTn(nc), Gen2(nc, nc)
      real(pr) :: Ge_i, Gen_i(nc), GeT_i, GeT2_i, GeTn_i(nc), Gen2_i(nc, nc)
      real(pr) :: ln_gammas(nc)
      real(pr) :: Ge_aux, Ge_aux2

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

      ! setup UNIFAC model
      model = setup_unifac(molecules)

      ! Call all Ge and derivatives
      call excess_gibbs(model, n, T, Ge, GeT, GeT2, Gen, GeTn, Gen2)
      call excess_gibbs(model, n + [0.0_pr, 0.001_pr, 0.0_pr], T, Ge_aux)
      call excess_gibbs(model, n - [0.0_pr, 0.001_pr, 0.0_pr], T, Ge_aux2)

      print *, Gen2(2,:)
      print *, (Ge_aux - 2 * Ge + Ge_aux2) / 0.001**2

      ! Call Ge and derivatives individually
      call excess_gibbs(model, n, T, Ge_i)
      call excess_gibbs(model, n, T, GeT=GeT_i)
      call excess_gibbs(model, n, T, GeT2=GeT2_i)
      call excess_gibbs(model, n, T, Gen=Gen_i)
      call excess_gibbs(model, n, T, GeTn=GeTn_i)
      call excess_gibbs(model, n, T, Gen2=Gen2_i)

      ! Call GeModel class method
      call model%ln_activity_coefficient(n, T, ln_gammas)
      
      ! ========================================================================
      ! Test against Caleb Bell's implementation
      ! ------------------------------------------------------------------------
      ! Ge
      call check(error, abs(Ge / n_t - (-3.223992676822129_pr)) <= 1e-10)

      ! Gen
      call check(error, allclose(Gen, [10.53032277_pr,  -2.37758326_pr, -36.65748951_pr], 1e-8_pr))

      ! ln_gammas
      call check(error, allclose(Gen / R / T, [0.84433781_pr, -0.19063836_pr, -2.93925506_pr], 1e-8_pr))
      call check(error, allclose(ln_gammas, [0.84433781_pr, -0.19063836_pr, -2.93925506_pr], 1e-8_pr))

      ! Gen2
      call check(error, allclose(Gen2(1,:), [-0.75249927_pr,  0.13440904_pr,  0.56413529_pr] * R*T, 1e-6_pr))
      call check(error, allclose(Gen2(2,:), [ 0.13440904_pr,  0.34708386_pr, -2.69840507_pr] * R*T, 1e-6_pr))
      call check(error, allclose(Gen2(3,:), [ 0.56413529_pr, -2.69840507_pr, 17.76056492_pr] * R*T, 1e-6_pr))

      ! dln_gammas_dn
      call check(error, allclose(Gen2(1,:) /R/T, [-0.752499273_pr, 0.134409037_pr, 0.564135287_pr], 1e-6_pr))
      call check(error, allclose(Gen2(2,:) /R/T, [0.1344090359_pr, 0.3470838559_pr, -2.698405064_pr], 1e-6_pr))
      call check(error, allclose(Gen2(3,:) /R/T, [0.5641352889_pr, -2.698405071_pr, 17.760564919_pr], 1e-6_pr))

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
      call excess_gibbs(model, n, T, Ge=Ge_i, GeT=GeT_i)
      call check(error, abs(Ge - Ge_i) <= 1e-10)
      call check(error, abs(GeT - GeT_i) <= 1e-10)

      call excess_gibbs(model, n, T, Ge=Ge_i, GeT2=GeT2_i)
      call check(error, abs(Ge - Ge_i) <= 1e-10)
      call check(error, abs(GeT2 - GeT2_i) <= 1e-10)

      call excess_gibbs(model, n, T, Ge=Ge_i, Gen=Gen_i)
      call check(error, abs(Ge - Ge_i) <= 1e-10)
      call check(error, allclose(Gen, Gen_i, 1e-10_pr))

      call excess_gibbs(model, n, T, Ge=Ge_i, GeTn=GeTn_i)
      call check(error, abs(Ge - Ge_i) <= 1e-10)
      call check(error, allclose(GeTn, GeTn_i, 1e-10_pr))

      call excess_gibbs(model, n, T, Ge=Ge_i, Gen2=Gen2_i)
      call check(error, abs(Ge - Ge_i) <= 1e-10)
      call check(error, allclose(Gen2(1,:), Gen2_i(1,:), 1e-10_pr))
      call check(error, allclose(Gen2(2,:), Gen2_i(2,:), 1e-10_pr))
      call check(error, allclose(Gen2(3,:), Gen2_i(3,:), 1e-10_pr))

      ! Ge_T
      call excess_gibbs(model, n, T, GeT=GeT_i, GeT2=GeT2_i)
      call check(error, abs(GeT - GeT_i) <= 1e-10)
      call check(error, abs(GeT2 - GeT2_i) <= 1e-10)

      call excess_gibbs(model, n, T, GeT=GeT_i, Gen=Gen_i)
      call check(error, abs(GeT - GeT_i) <= 1e-10)
      call check(error, allclose(Gen, Gen_i, 1e-10_pr))

      call excess_gibbs(model, n, T, GeT=GeT_i, GeTn=GeTn_i)
      call check(error, abs(GeT - GeT_i) <= 1e-10)
      call check(error, allclose(GeTn, GeTn_i, 1e-10_pr))

      call excess_gibbs(model, n, T, GeT=GeT_i, Gen2=Gen2_i)
      call check(error, abs(GeT - GeT_i) <= 1e-10)
      call check(error, allclose(Gen2(1,:), Gen2_i(1,:), 1e-10_pr))
      call check(error, allclose(Gen2(2,:), Gen2_i(2,:), 1e-10_pr))
      call check(error, allclose(Gen2(3,:), Gen2_i(3,:), 1e-10_pr))

      ! Ge_T2
      call excess_gibbs(model, n, T, GeT2=GeT2_i, Gen=Gen_i)
      call check(error, abs(GeT2 - GeT2_i) <= 1e-10)
      call check(error, allclose(Gen, Gen_i, 1e-10_pr))

      call excess_gibbs(model, n, T, GeT2=GeT2_i, GeTn=GeTn_i)
      call check(error, abs(GeT2 - GeT2_i) <= 1e-10)
      call check(error, allclose(GeTn, GeTn_i, 1e-10_pr))

      call excess_gibbs(model, n, T, GeT2=GeT2_i, Gen2=Gen2_i)
      call check(error, abs(GeT2 - GeT2_i) <= 1e-10)
      call check(error, allclose(Gen2(1,:), Gen2_i(1,:), 1e-10_pr))
      call check(error, allclose(Gen2(2,:), Gen2_i(2,:), 1e-10_pr))
      call check(error, allclose(Gen2(3,:), Gen2_i(3,:), 1e-10_pr))

      ! Gen_i
      call excess_gibbs(model, n, T, Gen=Gen_i, GeTn=GeTn_i)
      call check(error, allclose(Gen, Gen_i, 1e-10_pr))
      call check(error, allclose(GeTn, GeTn_i, 1e-10_pr))

      call excess_gibbs(model, n, T, Gen=Gen_i, Gen2=Gen2_i)
      call check(error, allclose(Gen, Gen_i, 1e-10_pr))
      call check(error, allclose(Gen2(1,:), Gen2_i(1,:), 1e-10_pr))
      call check(error, allclose(Gen2(2,:), Gen2_i(2,:), 1e-10_pr))
      call check(error, allclose(Gen2(3,:), Gen2_i(3,:), 1e-10_pr))

      ! ========================================================================
      ! Just one triplet call test
      ! ------------------------------------------------------------------------
      call excess_gibbs(model, n, T, Ge=Ge_i, GeT=GeT_i, Gen2=Gen2_i)
      call check(error, abs(Ge - Ge_i) <= 1e-10)
      call check(error, abs(GeT - GeT_i) <= 1e-10)
      call check(error, allclose(Gen2(1,:), Gen2_i(1,:), 1e-10_pr))
      call check(error, allclose(Gen2(2,:), Gen2_i(2,:), 1e-10_pr))
      call check(error, allclose(Gen2(3,:), Gen2_i(3,:), 1e-10_pr))
   end subroutine test_against_caleb_thermo
end module test_unifac
