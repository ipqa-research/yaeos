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
         new_unittest("Bubble pressure", test_bubble_pressure), &
         new_unittest("Dew pressure", test_dew_pressure) &
         ]
   end subroutine collect_suite

   subroutine test_bubble_pressure(error)
      use yaeos, only: pr, EquilibriaState, saturation_pressure, ArModel
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

      bubble = saturation_pressure(model, n, T, kind="bubble")
      call check(error, maxval(abs(bubble%x - x)) < abs_tolerance)
      call check(error, maxval(abs(bubble%y - y)) < abs_tolerance)
      call check(error, abs(bubble%p - p) < abs_tolerance)
   end subroutine
   
   subroutine test_dew_pressure(error)
      use yaeos, only: pr, EquilibriaState, saturation_pressure, ArModel
      use yaeos__phase_equilibria_auxiliar, only: k_wilson
      use fixtures_models, only: binary_PR76
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 2

      class(ArModel), allocatable :: model
      type(EquilibriaState) :: dew

      real(pr) :: x(nc)  = [0.4, 0.6]
      real(pr) :: y(nc) = [0.84203837140677695, 0.15796162859550475]
      real(pr) :: P = 12.124711835829961     

      real(pr) :: n(nc), k(nc), t

      integer :: i

      n = [0.4_pr, 0.6_pr]
      T = 300
      model = binary_PR76()

      k = k_wilson(model, T, P)
      y = n * k
      do i=300,700
         t = real(i, pr)
         dew = saturation_pressure(model, n, T, kind="dew", p0=100._pr, y0=y)
         if (.not. any(isnan(dew%y))) y = dew%y
         print *, dew%t, dew%p, dew%iters
      end do
   end subroutine

end module test_saturation
