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
         new_unittest("Pure Psat", test_pure_psat), &
         new_unittest("Bubble pressure", test_bubble_pressure), &
         new_unittest("Dew pressure", test_dew_pressure), &
         new_unittest("Bubble temperature", test_bubble_temperature), &
         new_unittest("Dew temperature", test_dew_temperature), &
         new_unittest("PT2 envelope", test_pt2_envelope), &
         new_unittest("PX2 envelope", test_px2_envelope) &
         ]
   end subroutine collect_suite

   subroutine test_bubble_pressure(error)
      use yaeos, only: pr, EquilibriumState, saturation_pressure, ArModel
      use fixtures_models, only: binary_PR76
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 2

      class(ArModel), allocatable :: model
      type(EquilibriumState) :: bubble

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
   end subroutine test_bubble_pressure

   subroutine test_dew_pressure(error)
      use yaeos, only: pr, EquilibriumState, saturation_pressure, ArModel
      use yaeos, only: k_wilson
      use fixtures_models, only: binary_PR76
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 2

      class(ArModel), allocatable :: model
      type(EquilibriumState) :: dew

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
   end subroutine test_dew_pressure

   subroutine test_dew_temperature(error)
      use yaeos, only: pr, EquilibriumState, saturation_temperature, ArModel
      use yaeos, only: k_wilson
      use fixtures_models, only: binary_PR76
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 2

      class(ArModel), allocatable :: model
      type(EquilibriumState) :: dew

      real(pr) :: x(nc) = [6.7245630132141868E-002, 0.93275436999337613]
      real(pr) :: y(nc)  = [0.4, 0.6]
      real(pr) :: P = 10.867413040635611

      real(pr) :: n(nc), k(nc), t

      integer :: i

      n = [0.4_pr, 0.6_pr]
      T = 240
      model = binary_PR76()

      dew = saturation_temperature(model, n, P, kind="dew", t0=250._pr)

      print *, dew

      call check(error, abs(dew%P-P) < abs_tolerance)
      call check(error, abs(dew%T-T) < abs_tolerance)
      call check(error, maxval(abs(dew%x-x)) < abs_tolerance)
      call check(error, maxval(abs(dew%y-y)) < abs_tolerance)
   end subroutine test_dew_temperature

   subroutine test_bubble_temperature(error)
      use yaeos, only: pr, EquilibriumState, saturation_temperature, ArModel, saturation_pressure, PTEnvel2, pt_envelope_2ph
      use fixtures_models, only: binary_PR76
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 2

      class(ArModel), allocatable :: model
      type(EquilibriumState) :: bubble

      real(pr) :: x(nc)  = [0.4, 0.6]
      real(pr) :: y(nc) = [0.84203837140677695, 0.15796162859550475]
      real(pr) :: P = 12.124711835829961

      real(pr) :: n(nc), t

      n = [0.4_pr, 0.6_pr]
      T = 200
      model = binary_PR76()
      bubble = saturation_temperature(model, n, P, kind="bubble",t0=201._pr)
      call check(error, maxval(abs(bubble%x - x)) < abs_tolerance)
      call check(error, maxval(abs(bubble%y - y)) < abs_tolerance)
      call check(error, abs(bubble%p - p) < abs_tolerance)
      call check(error, abs(bubble%T - T) < abs_tolerance)
   end subroutine test_bubble_temperature

   subroutine test_pt2_envelope(error)
      use yaeos, only: pr, EquilibriumState, saturation_pressure, ArModel
      use yaeos, only: &
         pt_envelope_2ph, PTEnvel2
      use fixtures_models, only: binary_PR76
      type(error_type), allocatable, intent(out) :: error

      type(EquilibriumState) :: bubble
      class(ArModel), allocatable :: model
      type(PTEnvel2) :: envelope
      real(pr) :: z(2) = [0.4_pr, 0.6_pr]

      z = z/sum(z)

      model = binary_PR76()

      bubble = saturation_pressure(model, z, 200._pr, kind="bubble", p0=10._pr)
      envelope = pt_envelope_2ph(&
         model, z, bubble &
         )
      print *, envelope%cps(1)
      call check(error, size(envelope%cps) == 1)
   end subroutine test_pt2_envelope

   subroutine test_px2_envelope(error)
      use yaeos, only: pr, EquilibriumState, saturation_pressure, ArModel
      use yaeos, only: &
         px_envelope_2ph, PXEnvel2
      use fixtures_models, only: binary_PR76
      type(error_type), allocatable, intent(out) :: error

      type(EquilibriumState) :: bubble
      class(ArModel), allocatable :: model
      type(PXEnvel2) :: envelope
      real(pr) :: z(2) = [0.01_pr, 0.99_pr]
      real(pr) :: z_inj(2) = [1,  0]
      integer :: i

      z = z/sum(z)

      model = binary_PR76()

      bubble = saturation_pressure(model, z, T=270._pr, kind="bubble", p0=10._pr)
      envelope = px_envelope_2ph(&
         model, z0=z, first_point=bubble, alpha0=0.0_pr, z_injection=z_inj&
         )
      call check(error, size(envelope%cps) == 1)
   end subroutine test_px2_envelope

   subroutine test_pure_psat(error)
      use yaeos, only: pr, ArModel, Psat
      use fixtures_models, only: binary_PR76
      type(error_type), allocatable, intent(out) :: error
      class(ArModel), allocatable :: model

      integer :: i, j
      real(pr) :: T, Psats(2), Psats_val(2)

      T = 150
      model = binary_PR76()
      Psats_val = [260.37450286310201, 30.028551527997834]

         do i=1,2
            Psats(i) = Psat(model, i, T)
         end do
   ! call check(error, maxval(abs(Psats-Psats_val)) < abs_tolerance)
   end subroutine test_pure_psat
end module test_saturation
