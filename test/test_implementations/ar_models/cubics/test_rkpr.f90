module test_rkpr
   use yaeos, only: pr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use auxiliar_functions, only: allclose, rel_error
   implicit none

contains
   subroutine collect_suite(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("Test RKPR Ar consistency mix", test_rkpr_cons_mixture), &
         new_unittest("Test RKPR Ar consistency pure", test_rkpr_cons_pure) &
         ]
   end subroutine collect_suite

   subroutine test_rkpr_cons_mixture(error)
      use yaeos, only: pr, RKPR, ArModel
      use yaeos__consistency, only: numeric_ar_derivatives, ar_consistency
      type(error_type), allocatable, intent(out) :: error

      class(ArModel), allocatable :: model, model_kij
      real(pr) :: tc(4), pc(4), w(4), zc(4)

      real(pr) :: n(4), t, v

      real(pr) :: Ar, ArV, ArT, Arn(size(n)), ArV2, ArT2, ArTV
      real(pr) :: ArVn(size(n)), ArTn(size(n)), Arn2(size(n), size(n))

      real(pr) :: Ar_num, ArV_num, ArT_num, Arn_num(size(n)), ArV2_num, ArT2_num
      real(pr) :: ArTV_num, ArVn_num(size(n)), ArTn_num(size(n))
      real(pr) :: Arn2_num(size(n), size(n))
      real(pr) :: eq31, eq33(size(n), size(n)), eq34(size(n)), eq36, eq37

      real(pr) :: kij(size(n), size(n)), lij(size(n), size(n))

      n = [1.5, 0.2, 0.7, 2.3]
      tc = [369.83, 425.12, 507.6, 540.2]
      pc = [42.48, 37.96, 30.25, 27.4]
      w = [0.152291, 0.200164, 0.301261, 0.349469]
      zc = [0.276, 0.274, 0.266, 0.261]

      t = 600_pr
      v = 0.5_pr

      ! ========================================================================
      ! Model without kij and lij
      ! ------------------------------------------------------------------------
      model = RKPR(tc, pc, w, zc)

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
      call check(error, rel_error(Ar, Ar_num) < 1e-5)
      call check(error, rel_error(ArV, ArV_num) < 1e-5)
      call check(error, rel_error(ArT, ArT_num) < 1e-5)
      call check(error, allclose(Arn, Arn_num, 1e-5_pr))
      call check(error, rel_error(ArV2, ArV2_num) < 1e-5)
      call check(error, rel_error(ArT2, ArT2_num) < 1e-5)
      call check(error, rel_error(ArTV, ArTV_num) < 1e-5)
      call check(error, allclose(ArVn, ArVn_num, 1e-5_pr))
      call check(error, allclose(ArTn, ArTn_num, 1e-5_pr))
      call check(error, maxval(rel_error(Arn2, Arn2_num)) < 1e-5)

      ! Consistency tests
      call check(error, abs(eq31) <= 1e-13)
      call check(error, maxval(abs(eq33)) < 1e-13)
      call check(error, maxval(abs(eq34)) < 1e-13)
      call check(error, abs(eq36) <= 1e-13)
      call check(error, abs(eq37) <= 1e-13)

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


      model_kij = RKPR(tc, pc, w, zc, kij, lij)

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
      call check(error, rel_error(Ar, Ar_num) < 1e-5)
      call check(error, rel_error(ArV, ArV_num) < 1e-5)
      call check(error, rel_error(ArT, ArT_num) < 1e-5)
      call check(error, allclose(Arn, Arn_num, 1e-5_pr))
      call check(error, rel_error(ArV2, ArV2_num) < 1e-5)
      call check(error, rel_error(ArT2, ArT2_num) < 1e-5)
      call check(error, rel_error(ArTV, ArTV_num) < 1e-5)
      call check(error, allclose(ArVn, ArVn_num, 1e-5_pr))
      call check(error, allclose(ArTn, ArTn_num, 1e-5_pr))
      call check(error, maxval(rel_error(Arn2, Arn2_num)) < 1e-5)

      ! Consistency tests
      call check(error, abs(eq31) <= 1e-13)
      call check(error, maxval(abs(eq33)) < 1e-14)
      call check(error, maxval(abs(eq34)) < 1e-13)
      call check(error, abs(eq36) <= 1e-14)
      call check(error, abs(eq37) <= 1e-14)
   end subroutine test_rkpr_cons_mixture

   subroutine test_rkpr_cons_pure(error)
      use yaeos, only: pr, RKPR, ArModel
      use yaeos__consistency, only: numeric_ar_derivatives, ar_consistency
      type(error_type), allocatable, intent(out) :: error

      class(ArModel), allocatable :: model
      real(pr) :: tc(1), pc(1), w(1), zc(1)

      real(pr) :: n(1), t, v

      real(pr) :: Ar, ArV, ArT, Arn(size(n)), ArV2, ArT2, ArTV
      real(pr) :: ArVn(size(n)), ArTn(size(n)), Arn2(size(n), size(n))

      real(pr) :: Ar_num, ArV_num, ArT_num, Arn_num(size(n)), ArV2_num, ArT2_num
      real(pr) :: ArTV_num, ArVn_num(size(n)), ArTn_num(size(n))
      real(pr) :: Arn2_num(size(n), size(n))
      real(pr) :: eq31, eq33(size(n), size(n)), eq34(size(n)), eq36, eq37

      n = [5.0]
      tc = [369.83]
      pc = [42.48]
      w = [0.152291]
      zc = [0.276]

      model = RKPR(tc, pc, w, zc)

      t = 600_pr
      v = 0.5_pr

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
      call check(error, rel_error(ArT2, ArT2_num) < 1e-5)
      call check(error, rel_error(ArTV, ArTV_num) < 1e-6)
      call check(error, allclose(ArVn, ArVn_num, 1e-6_pr))
      call check(error, allclose(ArTn, ArTn_num, 1e-6_pr))
      call check(error, maxval(rel_error(Arn2, Arn2_num)) < 1e-5)

      ! Consistency tests
      call check(error, abs(eq31) <= 1e-13)
      call check(error, maxval(abs(eq33)) < 1e-14)
      call check(error, maxval(abs(eq34)) < 1e-14)
      call check(error, abs(eq36) <= 1e-14)
      call check(error, abs(eq37) <= 1e-14)
   end subroutine test_rkpr_cons_pure

end module test_rkpr

