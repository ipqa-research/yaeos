program test_px3_solve_point_beta
   !! Three-phase Px point solver started from a poor phase-fraction guess.
   !!
   !! A three-phase point is first converged from the usual initialization
   !! (unstable bubble point + flash). Restarting the solver from the same
   !! compositions but with \(\beta\) a thousand times smaller forces the
   !! step limiter on \(\beta\); the solver must still converge to the same
   !! point, and the point must satisfy isofugacity between the three phases
   !! and the material balance. The same is checked for a restart from a
   !! shifted \(\alpha\) with the pressure specified (damped \(\alpha\) steps).
   use yaeos, only: pr, CubicEoS, EquilibriumState, saturation_pressure, flash
   use yaeos__equilibria_boundaries_phase_envelopes_px3, only: solve_point, get_values_from_X
   use fixtures_models, only: multicomponent_PR
   use testing_aux, only: assert, test_title
   implicit none

   integer, parameter :: nc = 12
   type(CubicEoS) :: eos
   type(EquilibriumState) :: sat
   type(EquilibriumState) :: fl
   real(pr) :: z0(nc)
   real(pr) :: zi(nc)
   real(pr) :: z(nc)
   real(pr) :: x(nc)
   real(pr) :: y(nc)
   real(pr) :: w(nc)
   real(pr) :: lnphi_x(nc)
   real(pr) :: lnphi_y(nc)
   real(pr) :: lnphi_w(nc)
   real(pr) :: a
   real(pr) :: P
   real(pr) :: beta
   real(pr) :: T
   real(pr) :: X_ref(2*nc + 3)
   real(pr) :: XX(2*nc + 3)
   real(pr) :: F(2*nc + 3)
   real(pr) :: dF(2*nc + 3, 2*nc + 3)
   integer :: ns
   integer :: its

   write(*, *) test_title("Three-phase Px point solver with a poor beta guess")

   eos = multicomponent_PR(z0, zi)
   a = 0.3_pr
   z = a*zi + (1 - a)*z0

   sat = saturation_pressure(eos, z, T=200._pr, kind="bubble", p0=20._pr)
   fl = flash(eos, z=z, T=sat%T + 50, P_spec=sat%P, iters=its)
   T = sat%T

   ! X = [ln Kx, ln Ky, ln P, alpha, beta], specify alpha
   XX = [log(fl%x/sat%y), log(fl%y/sat%y), log(sat%P), a, fl%beta]
   ns = 2*nc + 2
   call solve_point(eos, z0, zi, T, ns, XX(ns), XX, F, dF, its, 1000)
   call assert(its < 1000 .and. maxval(abs(F)) < 1e-5_pr, "Reference three-phase point converged")
   X_ref = XX

   ! Restart with a very small beta
   XX(2*nc + 3) = X_ref(2*nc + 3)*1e-3_pr
   call solve_point(eos, z0, zi, T, ns, XX(ns), XX, F, dF, its, 1000)
   call assert(its < 1000, "Poor-beta restart converged")
   call assert(its > 3, "Poor-beta restart needed the damped beta steps")
   call assert(abs(XX(2*nc + 3) - X_ref(2*nc + 3)) < 1e-4_pr, "Poor-beta restart: same beta")
   call assert(abs(XX(2*nc + 1) - X_ref(2*nc + 1)) < 1e-5_pr, "Poor-beta restart: same pressure")
   call assert(abs(XX(2*nc + 2) - a) < 1e-12_pr, "Poor-beta restart: alpha kept at its specification")

   ! Specify the pressure instead and start from a shifted alpha: the
   ! damped alpha steps must bring the point back to the reference one.
   XX = X_ref
   XX(2*nc + 2) = X_ref(2*nc + 2) + 0.15_pr
   ns = 2*nc + 1
   call solve_point(eos, z0, zi, T, ns, X_ref(ns), XX, F, dF, its, 1000)
   call assert(its < 1000, "Shifted-alpha restart (P specified) converged")
   call assert(abs(XX(2*nc + 2) - X_ref(2*nc + 2)) < 1e-4_pr, "Shifted-alpha restart: back to the reference alpha")
   call assert(abs(XX(2*nc + 3) - X_ref(2*nc + 3)) < 1e-4_pr, "Shifted-alpha restart: same beta")

   call get_values_from_X(z0, zi, XX, x, y, w, P, a, beta)

   ! Isofugacity of the three phases (stable roots, as in the solver)
   call eos%lnphi_pt(x, P, T, root_type="stable", lnPhi=lnphi_x)
   call eos%lnphi_pt(y, P, T, root_type="stable", lnPhi=lnphi_y)
   call eos%lnphi_pt(w, P, T, root_type="stable", lnPhi=lnphi_w)
   call assert(maxval(abs(log(x) + lnphi_x - log(w) - lnphi_w)) < 1e-4_pr, "Isofugacity between x and w")
   call assert(maxval(abs(log(y) + lnphi_y - log(w) - lnphi_w)) < 1e-4_pr, "Isofugacity between y and w")

   ! Material balance z = beta*y + (1-beta)*x
   z = a*zi + (1 - a)*z0
   call assert(maxval(abs(beta*y + (1 - beta)*x - z)) < 1e-6_pr, "Material balance of the two main phases")
   call assert(abs(sum(w) - 1) < 1e-5_pr, "Incipient phase is normalized")
end program test_px3_solve_point_beta
