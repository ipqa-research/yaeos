module test_cubic_alphas
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none

contains
   subroutine collect_suite(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("AlphaSoave", test_alpha_soave) &
                  ]
   end subroutine

   subroutine test_alpha_soave(error)
      use yaeos_constants, only: pr
      use yaeos_models_ar_cubic_alphas, only: AlphaSoave
      type(error_type), allocatable, intent(out) :: error

      type(AlphaSoave) :: alpha

      integer, parameter :: n = 2
      real(pr) :: Tr(n), k(n)
      real(pr) :: a(n), dadt(n), dadt2(n)

      real(pr) :: aval(n) = [1.1524213446238356, 1.0594365081389596]
      real(pr) :: davaldt(n) = [-0.33947331922020552,-0.14556349186104045]
      real(pr) :: davaldt2(n) = [0.47434164902525683,0.15556349186104046]     

      Tr = [0.4_pr, 0.5_pr]
      k = [0.2_pr, 0.1_pr]

      alpha%k = k
      call alpha%alpha(Tr, a, dadt, dadt2)

      call check(error, maxval(abs(a - aval)) < 1e-5)
      call check(error, maxval(abs(dadt - davaldt)) < 1e-5)
      call check(error, maxval(abs(dadt2 - davaldt2)) < 1e-5)
   end subroutine
end module
