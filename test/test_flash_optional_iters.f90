program test_flash_optional_iters
   !! Regression: `flash` and `flash_no_beta_limits` must work when the
   !! optional `iters` argument is omitted, and report the same iteration
   !! count in the result as when it is present.
   use yaeos, only: pr, ArModel, PengRobinson76, EquilibriumState, flash
   use yaeos__equilibria_flash, only: flash_no_beta_limits
   use testing_aux, only: assert, test_title
   implicit none
   integer, parameter :: nc = 2
   class(ArModel), allocatable :: model
   type(EquilibriumState) :: with_iters
   type(EquilibriumState) :: without_iters
   real(pr) :: z(nc)
   integer :: iters

   write(*, *) test_title("Flash with optional iters omitted")

   model = PengRobinson76( &
      [190.564_pr, 425.12_pr], [45.99_pr, 37.96_pr], [0.0115478_pr, 0.200164_pr])
   z = [0.3_pr, 0.7_pr]

   with_iters = flash(model, z, t=300.0_pr, p_spec=40.0_pr, iters=iters)
   without_iters = flash(model, z, t=300.0_pr, p_spec=40.0_pr)
   call assert(with_iters%kind == "split", "flash with iters finds a split")
   call assert(without_iters%kind == "split", "flash without iters finds a split")
   call assert(without_iters%iters == iters, "iteration count stored in the result")
   call assert(abs(without_iters%beta - with_iters%beta) < 1e-12_pr, "same beta with or without iters")

   with_iters = flash_no_beta_limits(model, z, t=300.0_pr, beta0=0.5_pr, p_spec=40.0_pr, iters=iters)
   without_iters = flash_no_beta_limits(model, z, t=300.0_pr, beta0=0.5_pr, p_spec=40.0_pr)
   call assert(without_iters%iters == iters, "flash_no_beta_limits: iteration count stored in the result")
end program test_flash_optional_iters
