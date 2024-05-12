module test_srk
   use yaeos, only: pr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use auxiliar_functions, only: allclose, rel_error
   implicit none

contains
   subroutine collect_suite(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("Test SRK Ar consistency mix", test_srk_cons_mixture), &
         new_unittest("Test SRK Ar consistency pure", test_srk_cons_pure) &
         ]
   end subroutine collect_suite

   subroutine test_srk_cons_mixture(error)
      use yaeos, only: pr, SoaveRedlichKwong, ArModel
      use yaeos_consistency, only: numeric_ar_derivatives, ar_consistency
      type(error_type), allocatable, intent(out) :: error

      class(ArModel), allocatable :: model, model_kij
      real(pr) :: tc(4), pc(4), w(4)

      real(pr) :: n(4), t, v

      real(pr) :: Ar, ArV, ArT, Arn(size(n)), ArV2, ArT2, ArTV
      real(pr) :: ArVn(size(n)), ArTn(size(n)), Arn2(size(n), size(n))

      real(pr) :: Ar_num, ArV_num, ArT_num, Arn_num(size(n)), ArV2_num, ArT2_num
      real(pr) :: ArTV_num, ArVn_num(size(n)), ArTn_num(size(n))
      real(pr) :: Arn2_num(size(n), size(n))
      real(pr) :: eq31, eq33(size(n), size(n)), eq34(size(n)), eq36, eq37

      real(pr) :: kij(size(n), size(n)), lij(size(n), size(n))

      n = [1.5, 0.2, 0.7, 2.3]
      tc = [190.564, 425.12, 300.11, 320.25]
      pc = [45.99, 37.96, 39.23, 40.21]
      w = [0.0115478, 0.200164, 0.3624, 0.298]

      t = 600_pr
      v = 0.5_pr

      ! ========================================================================
      ! Model without kij and lij
      ! ------------------------------------------------------------------------
      model = SoaveRedlichKwong(tc, pc, w)

      call model%residual_helmholtz(&
         n, v, t, Ar=Ar, ArV=ArV, ArT=ArT, Arn=Arn, &
         ArTV=ArTV, ArV2=ArV2, ArT2=ArT2, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2)

      call numeric_ar_derivatives(&
         model, n, v, t, d_n = 0.0001_pr, d_v = 0.0001_pr, d_t = 0.01_pr, &
         Ar=Ar_num, ArV=ArV_num, ArT=ArT_num, ArTV=ArTV_num, ArV2=ArV2_num, &
         ArT2=ArT2_num, Arn=Arn_num, ArVn=ArVn_num, ArTn=ArTn_num, &
         Arn2=Arn2_num &
         )

      call ar_consistency(&
         model, n, v, t, eq31=eq31, eq33=eq33, eq34=eq34, eq36=eq36, eq37=eq37 &
         )

      ! Numeric derivatives
      call check(error, rel_error(Ar, Ar_num) < 1e-6)
      call check(error, rel_error(ArV, ArV_num) < 1e-6)
      call check(error, rel_error(ArT, ArT_num) < 1e-6)
      call check(error, allclose(Arn, Arn_num, 1e-6_pr))
      call check(error, rel_error(ArV2, ArV2_num) < 1e-6)
      call check(error, rel_error(ArT2, ArT2_num) < 1e-6)
      call check(error, rel_error(ArTV, ArTV_num) < 1e-6)
      call check(error, allclose(ArVn, ArVn_num, 1e-6_pr))
      call check(error, allclose(ArTn, ArTn_num, 1e-6_pr))
      call check(error, maxval(rel_error(Arn2, Arn2_num)) < 1e-6)

      ! Consistency tests
      call check(error, abs(eq31) <= 1e-15)
      call check(error, maxval(abs(eq33)) < 1e-15)
      call check(error, maxval(abs(eq34)) < 1e-15)
      call check(error, abs(eq36) <= 1e-15)
      call check(error, abs(eq37) <= 1e-15)

      ! ========================================================================
      ! Model with kij and lij
      ! ------------------------------------------------------------------------
      kij = reshape([&
         0.0_pr, 0.1_pr, 0.2_pr, 0.1_pr, &
         0.1_pr, 0.0_pr, 0.3_pr, 0.25_pr, &
         0.2_pr, 0.3_pr, 0.0_pr, 0.18_pr, &
         0.1_pr, 0.25_pr, 0.18_pr, 0.0_pr], [size(n), size(n)])

      lij = reshape([&
         0.0_pr, 0.001_pr, 0.002_pr, 0.001_pr, &
         0.001_pr, 0.0_pr, 0.003_pr, 0.0025_pr, &
         0.002_pr, 0.003_pr, 0.0_pr, 0.0018_pr, &
         0.001_pr, 0.0025_pr, 0.0018_pr, 0.0_pr], [size(n), size(n)])


      model_kij = SoaveRedlichKwong(tc, pc, w, kij, lij)

      call model_kij%residual_helmholtz(&
         n, v, t, Ar=Ar, ArV=ArV, ArT=ArT, Arn=Arn, &
         ArTV=ArTV, ArV2=ArV2, ArT2=ArT2, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2)

      call numeric_ar_derivatives(&
         model_kij, n, v, t, d_n = 0.0001_pr, d_v = 0.0001_pr, d_t = 0.01_pr, &
         Ar=Ar_num, ArV=ArV_num, ArT=ArT_num, ArTV=ArTV_num, ArV2=ArV2_num, &
         ArT2=ArT2_num, Arn=Arn_num, ArVn=ArVn_num, ArTn=ArTn_num, &
         Arn2=Arn2_num &
         )

      call ar_consistency(&
         model_kij, n, v, t, eq31=eq31, eq33=eq33, &
         eq34=eq34, eq36=eq36, eq37=eq37 &
         )

      ! Numeric derivatives
      call check(error, rel_error(Ar, Ar_num) < 1e-6)
      call check(error, rel_error(ArV, ArV_num) < 1e-6)
      call check(error, rel_error(ArT, ArT_num) < 1e-6)
      call check(error, allclose(Arn, Arn_num, 1e-6_pr))
      call check(error, rel_error(ArV2, ArV2_num) < 1e-6)
      call check(error, rel_error(ArT2, ArT2_num) < 1e-6)
      call check(error, rel_error(ArTV, ArTV_num) < 1e-6)
      call check(error, allclose(ArVn, ArVn_num, 1e-6_pr))
      call check(error, allclose(ArTn, ArTn_num, 1e-6_pr))
      call check(error, maxval(rel_error(Arn2, Arn2_num)) < 1e-6)

      ! Consistency tests
      call check(error, abs(eq31) <= 1e-14)
      call check(error, maxval(abs(eq33)) < 1e-15)
      call check(error, maxval(abs(eq34)) < 1e-15)
      call check(error, abs(eq36) <= 1e-15)
      call check(error, abs(eq37) <= 1e-15)
   end subroutine test_srk_cons_mixture

   subroutine test_srk_cons_pure(error)
      use yaeos, only: pr, SoaveRedlichKwong, ArModel
      use yaeos_consistency, only: numeric_ar_derivatives, ar_consistency
      type(error_type), allocatable, intent(out) :: error

      class(ArModel), allocatable :: model
      real(pr) :: tc(1), pc(1), w(1)

      real(pr) :: n(1), t, v

      real(pr) :: Ar, ArV, ArT, Arn(size(n)), ArV2, ArT2, ArTV
      real(pr) :: ArVn(size(n)), ArTn(size(n)), Arn2(size(n), size(n))

      real(pr) :: Ar_num, ArV_num, ArT_num, Arn_num(size(n)), ArV2_num, ArT2_num
      real(pr) :: ArTV_num, ArVn_num(size(n)), ArTn_num(size(n))
      real(pr) :: Arn2_num(size(n), size(n))

      real(pr) :: eq31, eq33(size(n), size(n)), eq34(size(n)), eq36, eq37

      n = [5.0]
      tc = [190.564]
      pc = [45.99]
      w = [0.0115478]

      t = 600_pr
      v = 0.5_pr

      model = SoaveRedlichKwong(tc, pc, w)

      call model%residual_helmholtz(&
         n, v, t, Ar=Ar, ArV=ArV, ArT=ArT, Arn=Arn, &
         ArTV=ArTV, ArV2=ArV2, ArT2=ArT2, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2)

      call numeric_ar_derivatives(&
         model, n, v, t, d_n = 0.0001_pr, d_v = 0.0001_pr, d_t = 0.01_pr, &
         Ar=Ar_num, ArV=ArV_num, ArT=ArT_num, ArTV=ArTV_num, ArV2=ArV2_num, &
         ArT2=ArT2_num, Arn=Arn_num, ArVn=ArVn_num, ArTn=ArTn_num, &
         Arn2=Arn2_num &
         )

      call ar_consistency(&
         model, n, v, t, eq31=eq31, eq33=eq33, eq34=eq34, eq36=eq36, eq37=eq37 &
         )

      ! Numeric derivatives
      call check(error, rel_error(Ar, Ar_num) < 1e-6)
      call check(error, rel_error(ArV, ArV_num) < 1e-6)
      call check(error, rel_error(ArT, ArT_num) < 1e-6)
      call check(error, allclose(Arn, Arn_num, 1e-6_pr))
      call check(error, rel_error(ArV2, ArV2_num) < 1e-6)
      call check(error, rel_error(ArT2, ArT2_num) < 1e-6)
      call check(error, rel_error(ArTV, ArTV_num) < 1e-6)
      call check(error, allclose(ArVn, ArVn_num, 1e-6_pr))
      call check(error, allclose(ArTn, ArTn_num, 1e-6_pr))
      call check(error, maxval(rel_error(Arn2, Arn2_num)) < 1e-6)

      ! Consistency tests
      call check(error, abs(eq31) <= 1e-14)
      call check(error, maxval(abs(eq33)) < 1e-15)
      call check(error, maxval(abs(eq34)) < 1e-15)
      call check(error, abs(eq36) <= 1e-15)
      call check(error, abs(eq37) <= 1e-15)
   end subroutine test_srk_cons_pure

end module test_srk
