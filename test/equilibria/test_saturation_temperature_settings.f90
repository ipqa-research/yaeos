program test_saturation_temperature_settings
   !! Regression: `saturation_temperature` must not change the module setting
   !! `iters_first_step` that `saturation_pressure` uses.
   use yaeos, only: pr, ArModel, PengRobinson76, EquilibriumState, saturation_temperature
   use yaeos__equilibria_saturation_points, only: iters_first_step
   use testing_aux, only: assert, test_title
   implicit none
   class(ArModel), allocatable :: model
   type(EquilibriumState) :: sat
   integer :: before

   write(*, *) test_title("saturation_temperature leaves iters_first_step alone")
   model = PengRobinson76( &
      [190.564_pr, 425.12_pr], [45.99_pr, 37.96_pr], [0.0115478_pr, 0.200164_pr])
   before = iters_first_step
   sat = saturation_temperature(model, [0.3_pr, 0.7_pr], 20.0_pr, kind="bubble", t0=250.0_pr)
   call assert(sat%kind == "bubble", "bubble temperature converges")
   call assert(iters_first_step == before, "iters_first_step unchanged by saturation_temperature")
end program test_saturation_temperature_settings
