module test_math
   use yaeos__constants, only: pr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use auxiliar_functions, only: allclose
   implicit none

contains

   subroutine collect_suite(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("Test dx_to_dn", test_dx_to_dn), &
         new_unittest("Test sq_error", test_sq_error), &
         new_unittest("Test cardano", test_cardano_method) &
         ]
   end subroutine collect_suite

   subroutine test_dx_to_dn(error)
      use yaeos__math, only: dx_to_dn
      type(error_type), allocatable, intent(out) :: error

      real(pr) :: dns(3)

      dns = dx_to_dn(&
         [0.2_pr, 0.7_pr, 0.1_pr], &
         [1053.0323_pr, -237.758335_pr, -3665.7490_pr]&
         )

      call check(error, allclose(dns, [1375.4315_pr, 84.6409421_pr, -3343.3497_pr], 1e-4_pr))
   end subroutine test_dx_to_dn

   subroutine test_sq_error(error)
      use yaeos__math, only: sq_error

      type(error_type), allocatable, intent(out) :: error

      real(pr) :: sim(5), exps(5), errors_sq(5)

      sim = [0.5_pr, 1.5_pr, 0.4_pr, 1.85_pr, 1.66_pr]
      exps = [0.45_pr, 1.51_pr, 0.42_pr, 1.88_pr, 1.68_pr]

      errors_sq = sq_error(exps, sim)

      call check(error, allclose(errors_sq, (sim - exps)**2, 1e10_pr))

   end subroutine test_sq_error

   subroutine test_cardano_method(error)
      use yaeos__math_linalg, only: pr, cubic_roots
      type(error_type), allocatable, intent(out) :: error
      real(pr) :: p(4)
      real(pr) :: rr(3)
      complex(pr) :: cr(3)
      real(pr) :: num
      integer :: flag

      p = [1, -10, 35, -50]
      call cubic_roots(p, rr, cr, flag)
      call check(error, flag == 1)
      call check(error, abs(rr(1) - 5.0_pr) < 1e-10_pr)

      p = [0.1, -2.6, 1., 1.]
      call cubic_roots(p, rr, cr, flag)
      call check(error, flag == -1)
      call check(error, &
         maxval(abs(rr - [-0.454216, 0.8601986, 25.594016])) < 1e-5_pr &
      )
   end subroutine test_cardano_method
end module test_math


