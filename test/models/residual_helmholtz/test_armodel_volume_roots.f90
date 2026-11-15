program test_armodel_volume_roots
   !! Regression for the generic `ArModel%volume` solver (models that do not
   !! override it, here a hyperdual PR76 checked against the CubicEoS solver):
   !! - "stable" must return the converged root when the other one fails;
   !! - the vapor initial guess must use the total moles, so volumes scale
   !!   with n.
   use yaeos, only: pr, ArModel
   use fixtures_models, only: binary_PR76_hd, binary_PR76
   use testing_aux, only: assert, test_title
   implicit none
   class(ArModel), allocatable :: hd
   class(ArModel), allocatable :: cub
   real(pr) :: z(2)
   real(pr) :: V
   real(pr) :: V_ref

   write(*, *) test_title("Generic ArModel volume: stable root and total moles")
   hd = binary_PR76_hd()
   cub = binary_PR76()
   z = [0.3_pr, 0.7_pr]

   ! 200 K, 100 bar: the liquid root exists; the solve from the ideal-gas
   ! guess does not converge. "stable" must still return the liquid root.
   call cub%volume(z, 100.0_pr, 200.0_pr, V_ref, root_type="stable")
   call hd%volume(z, 100.0_pr, 200.0_pr, V, root_type="stable")
   call assert(V > 0.0_pr, "stable root found when only the liquid solve converges")
   call assert(abs(V - V_ref) < 1e-6_pr*V_ref, "stable root equals the CubicEoS stable root")

   ! Vapor root for 5 mol at 200 K and 1 bar: five times the molar volume.
   call cub%volume(z, 1.0_pr, 200.0_pr, V_ref, root_type="vapor")
   call hd%volume(5*z, 1.0_pr, 200.0_pr, V, root_type="vapor")
   call assert(abs(V - 5*V_ref) < 1e-6_pr*5*V_ref, "vapor volume scales with total moles")
end program test_armodel_volume_roots
