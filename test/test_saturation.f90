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
         new_unittest("Dew pressure", test_dew_pressure), &
         new_unittest("Bubble temperature", test_bubble_temperature), &
         new_unittest("Dew temperature", test_dew_temperature), &
         new_unittest("PT envelope", test_envelope) &
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

      real(pr) :: x(nc) = [6.7245630132141868E-002, 0.93275436999337613]
      real(pr) :: y(nc)  = [0.4, 0.6]
      real(pr) :: P = 10.867413040635611

      real(pr) :: n(nc), k(nc), t

      integer :: i

      n = [0.4_pr, 0.6_pr]
      T = 240
      model = binary_PR76()

      k = k_wilson(model, T, P)
      dew = saturation_pressure(model, n, T, kind="dew", p0=17._pr)
      call check(error, abs(dew%P-P) < abs_tolerance)
      call check(error, abs(dew%T-T) < abs_tolerance)
      call check(error, maxval(abs(dew%x-x)) < abs_tolerance)
      call check(error, maxval(abs(dew%y-y)) < abs_tolerance)
   end subroutine

   subroutine test_dew_temperature(error)
      use yaeos, only: pr, EquilibriaState, saturation_temperature, ArModel
      use yaeos__phase_equilibria_auxiliar, only: k_wilson
      use fixtures_models, only: binary_PR76
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 2

      class(ArModel), allocatable :: model
      type(EquilibriaState) :: dew

      real(pr) :: x(nc) = [6.7245630132141868E-002, 0.93275436999337613]
      real(pr) :: y(nc)  = [0.4, 0.6]
      real(pr) :: P = 10.867413040635611

      real(pr) :: n(nc), k(nc), t

      integer :: i

      n = [0.4_pr, 0.6_pr]
      T = 240
      model = binary_PR76()

      dew = saturation_temperature(model, n, P, kind="dew", t0=250._pr)
      call check(error, abs(dew%P-P) < abs_tolerance)
      call check(error, abs(dew%T-T) < abs_tolerance)
      call check(error, maxval(abs(dew%x-x)) < abs_tolerance)
      call check(error, maxval(abs(dew%y-y)) < abs_tolerance)
   end subroutine
   
   subroutine test_bubble_temperature(error)
      use yaeos, only: pr, EquilibriaState, saturation_temperature, ArModel
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

      bubble = saturation_temperature(model, n, P, kind="bubble", t0=100._pr)
      call check(error, maxval(abs(bubble%x - x)) < abs_tolerance)
      call check(error, maxval(abs(bubble%y - y)) < abs_tolerance)
      call check(error, abs(bubble%p - p) < abs_tolerance)
      call check(error, abs(bubble%T - T) < abs_tolerance)
   end subroutine

   subroutine test_envelope(error)
      use yaeos, only: pr, EquilibriaState, saturation_pressure, ArModel
      use yaeos__phase_equilibria_boundaries_phase_envelopes_pt, only: &
         pt_envelope_2ph, PTEnvel2
      use fixtures_models, only: binary_PR76
      type(error_type), allocatable, intent(out) :: error

      type(EquilibriaState) :: bubble
      class(ArModel), allocatable :: model
      type(PTEnvel2) :: envelope
      real(pr) :: z(2) = [0.4_pr, 0.6_pr]
      integer :: i

      z = z/sum(z)

      model = binary_PR76()

      bubble = saturation_pressure(model, z, 200._pr, kind="bubble", p0=10._pr)
      envelope = pt_envelope_2ph(&
         model, z, bubble &
      )
      call check(error, size(envelope%cps) == 1)
   end subroutine
end module test_saturation
