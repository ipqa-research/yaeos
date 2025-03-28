program main
   use yaeos
   use yaeos__equilibria_boundaries_phase_envelopes_pt3, only: &
      solve_point, PTEnvel3, pt_envelope_3ph, get_values_from_X
   use fixtures_models, only: multicomponent_PR, asphaltenes_srk
   use testing_aux, only: assert, test_title
   implicit none

   type(CubicEoS) :: eos

   type(PTEnvel2) :: env
   type(PTEnvel3) :: env3

   type(EquilibriumState) :: sat, fl
   integer, parameter :: nc=12
   real(pr) :: z0(nc), zi(nc), a
   real(pr) :: z(nc), P, T

   real(pr) :: w(nc), y(nc), x(nc), lnkx(nc), lnky(nc), lnP, lnT, beta
   real(pr) :: XX(2*nc+3), F(2*nc+3), dF(2*nc+3, 2*nc+3)

   integer :: ns, its, i

   write(*, *) test_title("3ph PT envelope")

   eos = multicomponent_PR(z0, zi)

   ! Define the composition
   a = 0.3
   z = a * zi + (1-a)*z0

   ! Calculation of unstable bubble line
   sat = saturation_pressure(eos, z, T=200._pr, kind="bubble", p0=20._pr)

   ! Calculation of a flash to initialize point
   fl = flash(eos, z=z, T=sat%T+50, P_spec=sat%P, iters=its)
   env = pt_envelope_2ph(eos, z, sat, points=1000)

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
   XX = [lnKx, lnKy, lnP, lnT, beta]

   ! Specify temperature
   ns = 2*nc + 2

   ! Solve point
   call solve_point(eos, z, ns, XX(ns), XX, F, dF, its, 1000)

   ! Obtain the values of the point
   call get_values_from_X(z, XX, x, y, w, P, T, beta)

   call assert(its < 1000, "solve_point convergence")
   call assert(abs(beta - 0.22) < 1e-2, "beta value")
   call assert(abs(P - 27.158429541826763) < 1e-2, "P value")
   call assert(abs(T - 200.) < 1e-2, "T value")

   ! Calculate the PT envelope using the converged point as the initial point
   env3 = pt_envelope_3ph(&
      eos, z=z, x0=x, y0=y, w0=w, beta0=beta, &
      P0=P, T0=T, ns0=ns, dS0=1e-5_pr, points=900)

   i = size(env3%P)

   print *, i, env3%P(i), env3%T(i), env3%beta(i)
   call assert(env3%T(i) < 220._pr, "Stop temperature of the PT envelope")
   call assert(env3%P(i) < 13._pr, "Stop pressure of the PT envelope")
   call assert(env3%beta(i) > 0.4_pr, "Stop beta of the PT envelope")
contains

   subroutine init_bub_with_asphaltenes(model, z, w, lnK, beta)
      !! Correction for fluids that contain asphaltenes, asuming that the
      !! asphaltenes are the last component of the composition vector.
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: z(:), w(:)
      real(pr), intent(out) :: lnK(:)
      real(pr), intent(out) :: beta
      real(pr) :: phi_w(size(z)), phi_y(size(z))

      y = 0
      y(nc) = 1
      call eos%lnphi_pt(y, sat%P, sat%T, root_type="liquid", lnphi=phi_y)
      call eos%lnphi_pt(w, sat%P, sat%T, root_type="vapor", lnphi=phi_w)
      lnK = phi_w - phi_y
      beta = z(nc)
   end subroutine init_bub_with_asphaltenes
end program main
