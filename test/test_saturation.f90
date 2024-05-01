module test_saturation
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use yaeos, only: pr
   implicit none

   real(pr) :: abs_tolerance = 1e-5_pr

contains
   subroutine collect_suite(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("Bubble pressure", test_bubble_pressure) &
         ]
   end subroutine collect_suite

   subroutine test_bubble_pressure(error)
      use yaeos, only: pr, EquilibriaState, bubble_pressure, ArModel
      use fixtures_models, only: binary_PR76
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 2

      class(ArModel), allocatable :: model
      type(EquilibriaState) :: bubble

      real(pr) :: x(nc)  = [0.4, 0.6]
      real(pr) :: y(nc) = [0.84203837140677695, 0.15796162859550475]
      real(pr) :: P = 12.124711835829961     

      real(pr) :: n(nc), t

      n = [0.4_pr, 0.6_pr]
      T = 200
      model = binary_PR76()

      bubble = bubble_pressure(model, n, T)
      call check(error, maxval(abs(bubble%x - x)) < abs_tolerance)
      call check(error, maxval(abs(bubble%y - y)) < abs_tolerance)
      call check(error, abs(bubble%p - p) < abs_tolerance)
   end subroutine

end module test_saturation
