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
         new_unittest("Test PSRK Ar consistency mix", test_psrk_cons_mixture), &
         new_unittest("Test PSRK Ar consistency pure", test_psrk_cons_pure) &
         ]
   end subroutine collect_suite

   subroutine test_psrk_cons_mixture(error)
      use yaeos, only: pr, PSRK, ArModel, Groups
      use yaeos__consistency, only: numeric_ar_derivatives, ar_consistency
      type(error_type), allocatable, intent(out) :: error

      class(ArModel), allocatable :: model
      real(pr) :: tc(2), pc(2), w(2), C(3, 2)
      type(Groups) :: molecules(2)

      real(pr) :: n(2), t, v

      real(pr) :: Ar, ArV, ArT, Arn(size(n)), ArV2, ArT2, ArTV
      real(pr) :: ArVn(size(n)), ArTn(size(n)), Arn2(size(n), size(n))

      real(pr) :: Ar_num, ArV_num, ArT_num, Arn_num(size(n)), ArV2_num, ArT2_num
      real(pr) :: ArTV_num, ArVn_num(size(n)), ArTn_num(size(n))
      real(pr) :: Arn2_num(size(n), size(n))
      real(pr) :: eq31, eq33(size(n), size(n)), eq34(size(n)), eq36, eq37

      n = [60, 40]
      tc = [304.21_pr, 553.8_pr]
      pc = [7.383e6_pr, 4.080358e6_pr] / 1e5
      w = [0.223621_pr, 0.213_pr]
      c(:, 1) = [0.8255_pr, 0.16755_pr, -1.7039_pr]
      c(:, 2)= [0.84082_pr, -0.39847_pr, 0.94148_pr]

      molecules(1)%groups_ids = [117]
      molecules(1)%number_of_groups = [2]
      molecules(2)%groups_ids = [2]
      molecules(2)%number_of_groups = [7]

      t = 600_pr
      v = 10.0_pr

      ! ========================================================================
      ! Model without kij and lij
      ! ------------------------------------------------------------------------
      model = PSRK(tc, pc, w, molecules, c1=C(1, :), c2=C(2, :), c3=C(3, :))

      call model%residual_helmholtz(&
         n, v, t, Ar=Ar, ArV=ArV, ArT=ArT, Arn=Arn, &
         ArTV=ArTV, ArV2=ArV2, ArT2=ArT2, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2)

      call numeric_ar_derivatives(&
         model, n, v, t, d_n = 0.01_pr, d_v = 0.01_pr, d_t = 0.01_pr, &
         Ar=Ar_num, ArV=ArV_num, ArT=ArT_num, ArTV=ArTV_num, ArV2=ArV2_num, &
         ArT2=ArT2_num, Arn=Arn_num, ArVn=ArVn_num, ArTn=ArTn_num, &
         Arn2=Arn2_num &
         )

      call ar_consistency(&
         model, n, v, t, eq31=eq31, eq33=eq33, eq34=eq34, eq36=eq36, eq37=eq37 &
         )

      ! Numeric derivatives
      call check(error, rel_error(Ar, Ar_num) < 1e-4)
      call check(error, rel_error(ArV, ArV_num) < 1e-3)
      call check(error, rel_error(ArT, ArT_num) < 1e-4)
      call check(error, allclose(Arn, Arn_num, 1e-4_pr))
      call check(error, rel_error(ArV2, ArV2_num) < 1e-3)
      call check(error, rel_error(ArT2, ArT2_num) < 1e-4)
      call check(error, rel_error(ArTV, ArTV_num) < 1e-3)
      call check(error, allclose(ArVn, ArVn_num, 1e-4_pr))
      call check(error, allclose(ArTn, ArTn_num, 1e-4_pr))
      call check(error, maxval(rel_error(Arn2, Arn2_num)) < 1e-4)

      ! Consistency tests
      call check(error, abs(eq31) <= 1e-13)
      call check(error, maxval(abs(eq33)) < 1e-13)
      call check(error, maxval(abs(eq34)) < 1e-13)
      call check(error, abs(eq36) <= 1e-13)
      call check(error, abs(eq37) <= 1e-13)
   end subroutine test_psrk_cons_mixture

   subroutine test_psrk_cons_pure(error)
      use yaeos, only: pr, PSRK, ArModel, Groups
      use yaeos__consistency, only: numeric_ar_derivatives, ar_consistency
      type(error_type), allocatable, intent(out) :: error

      class(ArModel), allocatable :: model
      type(Groups) :: molecules(1)
      real(pr) :: tc(1), pc(1), w(1), c(1, 3)

      real(pr) :: n(1), t, v

      real(pr) :: Ar, ArV, ArT, Arn(size(n)), ArV2, ArT2, ArTV
      real(pr) :: ArVn(size(n)), ArTn(size(n)), Arn2(size(n), size(n))

      real(pr) :: Ar_num, ArV_num, ArT_num, Arn_num(size(n)), ArV2_num, ArT2_num
      real(pr) :: ArTV_num, ArVn_num(size(n)), ArTn_num(size(n))
      real(pr) :: Arn2_num(size(n), size(n))

      real(pr) :: eq31, eq33(size(n), size(n)), eq34(size(n)), eq36, eq37

      n = [5]
      tc = [304.21_pr]
      pc = [7.383e6_pr] / 1e5
      w = [0.223621_pr]
      c(1, :) = [0.8255_pr, 0.16755_pr, -1.7039_pr]

      molecules(1)%groups_ids = [117]
      molecules(1)%number_of_groups = [2]

      t = 600_pr
      v = 15._pr

      model = PSRK(tc, pc, w, molecules, c1=c(:, 1), c2=c(:, 2), c3=c(:, 3))

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
      call check(error, rel_error(Ar, Ar_num) < 1e-4)
      call check(error, rel_error(ArV, ArV_num) < 1e-3)
      call check(error, rel_error(ArT, ArT_num) < 1e-4)
      call check(error, allclose(Arn, Arn_num, 1e-4_pr))
      call check(error, rel_error(ArV2, ArV2_num) < 1e-3)
      call check(error, rel_error(ArT2, ArT2_num) < 1e-4)
      call check(error, rel_error(ArTV, ArTV_num) < 1e-3)
      call check(error, allclose(ArVn, ArVn_num, 1e-4_pr))
      call check(error, allclose(ArTn, ArTn_num, 1e-4_pr))
      call check(error, maxval(rel_error(Arn2, Arn2_num)) < 1e-4)

      ! Consistency tests
      call check(error, abs(eq31) <= 1e-13)
      call check(error, maxval(abs(eq33)) < 1e-13)
      call check(error, maxval(abs(eq34)) < 1e-13)
      call check(error, abs(eq36) <= 1e-13)
      call check(error, abs(eq37) <= 1e-13)
   end subroutine test_psrk_cons_pure

end module test_psrk
