program main
   use yaeos
   use yaeos__equilibria_boundaries_phase_envelopes_px3, only: &
      solve_point, PXEnvel3, px_envelope_3ph, get_values_from_X
   use fixtures_models, only: multicomponent_PR, asphaltenes_srk
   use testing_aux, only: assert, test_title
   implicit none

   type(CubicEoS) :: eos

   type(PXEnvel3) :: env3

   type(EquilibriumState) :: sat, fl
   integer, parameter :: nc=12
   real(pr) :: z0(nc), zi(nc), a
   real(pr) :: z(nc), P, T

   real(pr) :: w(nc), y(nc), x(nc), lnkx(nc), lnky(nc), lnP, lnT, beta
   real(pr) :: XX(2*nc+3), F(2*nc+3), dF(2*nc+3, 2*nc+3)

   integer :: ns, its, i

   write(*, *) test_title("3ph PX envelope")

   eos = multicomponent_PR(z0, zi)

   ! Define the composition
   a = 0.3
   z = a * zi + (1-a)*z0

   ! Calculation of unstable bubble line
   sat = saturation_pressure(eos, z, T=200._pr, kind="bubble", p0=20._pr)

   ! Calculation of a flash to initialize point
   fl = flash(eos, z=z, T=sat%T+50, P_spec=sat%P, iters=its)

   ! The incipient phase is the unstable gas phase
   w = sat%y

   ! Initialize the main phases with the flash result
   x = fl%x
   y = fl%y

   ! Set up variables for the point solver
   lnP = log(sat%P)
   lnT = log(sat%T)
   lnkx = log(x/w)
   lnky = log(y/w)
   beta = fl%beta
   XX = [lnKx, lnKy, lnP, a, beta]

   ! Specify temperature
   ns = 2*nc + 2

   ! Solve point
   call solve_point(eos, z0, zi, sat%T, ns, XX(ns), XX, F, dF, its, 1000)

   ! Obtain the values of the point
   call get_values_from_X(z0, zi, XX, x, y, w, P, a, beta)

   call assert(its < 1000, "solve_point did not converge")
   call assert(abs(beta - 0.22) < 1e-2, "beta value")
   call assert(abs(P - 27.158429541826763) < 1e-2, "P value")
   call assert(abs(a - 0.3) < 1e-2, "alpha value")

   ! Calculate the PT envelope using the converged point as the initial point
   env3 = px_envelope_3ph(&
      eos, z0=z0, zi=zi, T=sat%T, x0=x, y0=y, w0=w, beta0=beta, &
      P0=P, a0=a, ns0=ns, dS0=0.01_pr, points=1000)

   i = size(env3%P)

   print *, "Number of points in the PT envelope: ", i
   print *, "Stop temperature of the PT envelope: ", env3%alpha(i)
   print *, "Stop pressure of the PT envelope: ", env3%P(i)
   print *, "Stop beta of the PT envelope: ", env3%beta(i)

   call assert(i > 10, "Number of points in the PT envelope")
   call assert(abs(env3%alpha(i) - 0.99) < 0.01, "Stop temperature of the PT envelope")
   call assert(abs(env3%P(i) - 2.58_pr) < 0.1, "Stop pressure of the PT envelope")
   call assert(abs(env3%beta(i) - 0.99) < 0.01_pr, "Stop beta of the PT envelope")

   call exit
end program main
