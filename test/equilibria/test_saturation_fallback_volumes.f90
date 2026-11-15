program test_saturation_fallback_volumes
   !! Regression: the full-Newton fallback of `saturation_pressure` and
   !! `saturation_temperature` must return molar volumes even when the input
   !! amounts do not sum to one (it used the raw amounts for the main phase).
   use yaeos, only: pr, ArModel, PengRobinson76, EquilibriumState, &
      saturation_pressure, saturation_temperature
   use yaeos__equilibria_saturation_points, only: iters_first_step
   use testing_aux, only: assert, test_title
   implicit none
   class(ArModel), allocatable :: model
   type(EquilibriumState) :: ref
   type(EquilibriumState) :: scaled
   real(pr) :: z(2)
   integer :: default_first_step

   write(*, *) test_title("Saturation fallback returns molar volumes")
   model = PengRobinson76( &
      [190.564_pr, 425.12_pr], [45.99_pr, 37.96_pr], [0.0115478_pr, 0.200164_pr])
   z = [0.3_pr, 0.7_pr]
   default_first_step = iters_first_step

   ! One scalar step only, so the full-Newton fallback produces the result.
   iters_first_step = 1
   ref = saturation_pressure(model, z, 300.0_pr, kind="bubble", p0=50.0_pr)
   scaled = saturation_pressure(model, 10*z, 300.0_pr, kind="bubble", p0=50.0_pr)
   iters_first_step = default_first_step
   call assert(ref%kind == "bubble" .and. scaled%kind == "bubble", "fallback bubble pressures converge")
   call assert(abs(scaled%P - ref%P) < 1e-6_pr*ref%P, "same pressure for n = z and n = 10 z")
   call assert(abs(scaled%Vx - ref%Vx) < 1e-8_pr*ref%Vx, "bubble pressure: Vx is molar for n = 10 z")
end program test_saturation_fallback_volumes
