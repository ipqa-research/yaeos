program test_saturation
   use yaeos, only: pr
   use testing_aux, only: assert, test_title
   implicit none

   real(pr) :: abs_tolerance = 1e-5_pr

   write(*, *) test_title("Saturation points and phase envelopes")
   call test_pure_psat()
   call test_bubble_pressure()
   call test_dew_pressure()
   call test_bubble_temperature()
   call test_dew_temperature()
   call test_pt2_envelope()
   call test_px2_envelope()
   write(*, *) ""

contains
   subroutine test_bubble_pressure()
      use yaeos, only: pr, EquilibriumState, saturation_pressure, ArModel
      use fixtures_models, only: binary_PR76

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
      call assert(maxval(abs(bubble%x - x)) < abs_tolerance, "Bubble pressure: x vs reference")
      call assert(maxval(abs(bubble%y - y)) < abs_tolerance, "Bubble pressure: y vs reference")
      call assert(abs(bubble%p - p) < abs_tolerance, "Bubble pressure: p vs reference")
   end subroutine test_bubble_pressure

   subroutine test_dew_pressure()
      use yaeos, only: pr, EquilibriumState, saturation_pressure, ArModel
      use yaeos, only: k_wilson
      use fixtures_models, only: binary_PR76

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
      call assert(abs(dew%P-P) < abs_tolerance, "Dew pressure: P vs reference")
      call assert(abs(dew%T-T) < abs_tolerance, "Dew pressure: T vs reference")
      call assert(maxval(abs(dew%x-x)) < abs_tolerance, "Dew pressure: x vs reference")
      call assert(maxval(abs(dew%y-y)) < abs_tolerance, "Dew pressure: y vs reference")
   end subroutine test_dew_pressure

   subroutine test_dew_temperature()
      use yaeos, only: pr, EquilibriumState, saturation_temperature, ArModel
      use yaeos, only: k_wilson
      use fixtures_models, only: binary_PR76

      integer, parameter :: nc = 2

      class(ArModel), allocatable :: model
      type(EquilibriumState) :: dew

      real(pr) :: x(nc) = [0.0673,  0.9327]
      real(pr) :: y(nc)  = [0.4, 0.6]
      real(pr) :: P = 10.867413040635611

      real(pr) :: n(nc), t

      n = [0.4_pr, 0.6_pr]
      T = 240
      model = binary_PR76()

      dew = saturation_temperature(model, n, P, kind="dew", t0=250._pr)

      call assert(abs(dew%P-P) < abs_tolerance, "Dew temperature: P vs reference")
      call assert(abs(dew%T-T) < abs_tolerance, "Dew temperature: T vs reference")
      call assert(maxval(abs(dew%x-x)) < 1e-4, "Dew temperature: x vs reference")
      call assert(maxval(abs(dew%y-y)) < 1e-4, "Dew temperature: y vs reference")
   end subroutine test_dew_temperature

   subroutine test_bubble_temperature()
      use yaeos, only: pr, EquilibriumState, saturation_temperature, ArModel, saturation_pressure, PTEnvel2, pt_envelope_2ph
      use fixtures_models, only: binary_PR76

      integer, parameter :: nc = 2

      class(ArModel), allocatable :: model
      type(EquilibriumState) :: bubble

      real(pr) :: x(nc)  = [0.4, 0.6]
      real(pr) :: y(nc) = [0.84230325312571430, 0.15819060946721245]
      real(pr) :: P = 12.124711835829961

      real(pr) :: n(nc), t

      n = [0.4_pr, 0.6_pr]
      T = 200
      model = binary_PR76()
      bubble = saturation_temperature(model, n, P, kind="bubble",t0=201._pr)
      call assert(maxval(abs(bubble%x - x)) < abs_tolerance, "Bubble temperature: x vs reference")
      call assert(maxval(abs(bubble%y - y)) < abs_tolerance, "Bubble temperature: y vs reference")
      call assert(abs(bubble%p - p) < abs_tolerance, "Bubble temperature: p vs reference")
      call assert(abs(bubble%T - T) < abs_tolerance, "Bubble temperature: T vs reference")
   end subroutine test_bubble_temperature

   subroutine test_pt2_envelope()
      use yaeos, only: pr, EquilibriumState, saturation_pressure, ArModel
      use yaeos, only: &
         pt_envelope_2ph, PTEnvel2
      use fixtures_models, only: binary_PR76

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
      call assert(size(envelope%cps) == 1, "PT2 envelope: exactly one critical point detected")
   end subroutine test_pt2_envelope

   subroutine test_px2_envelope()
      use yaeos, only: pr, EquilibriumState, saturation_pressure, ArModel
      use yaeos, only: &
         px_envelope_2ph, PXEnvel2
      use fixtures_models, only: binary_PR76

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
         model, z0=z, first_point=bubble, alpha0=0.0_pr, z_injection=z_inj, delta_0=0.001_pr&
         )

      call assert(size(envelope%cps) == 1, "PX2 envelope: exactly one critical point detected")
   end subroutine test_px2_envelope

   subroutine test_pure_psat()
      use yaeos, only: pr, ArModel
      use fixtures_models, only: binary_PR76
      class(ArModel), allocatable :: model

      integer :: i, j
      real(pr) :: T, Psats(2), Psats_val(2), Vv, Vl

      T = 150
      model = binary_PR76()
      Psats_val = [260.37450286310201, 30.028551527997834]

      do i=1,2
         Psats(i) = model%Psat_pure(i, T, Vl=Vl, Vv=Vv)
      end do
      ! call assert(maxval(abs(Psats-Psats_val)) < abs_tolerance, "Pure Psat: Psat vs reference")
   end subroutine test_pure_psat
end program test_saturation
