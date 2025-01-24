program main
   use yaeos
   use yaeos__equilibria_boundaries_phase_envelopes_pt3, only: solve_point, PTEnvel3, pt_envelope_3ph
   use fixtures_models, only: multicomponent_PR, asphaltenes_srk
   use testing_aux, only: assert, test_title
   implicit none

   type(CubicEoS) :: eos

   type(PTEnvel2) :: env
   type(PTEnvel3) :: env3

   type(EquilibriumState) :: sat, cp, fl
   integer, parameter :: nc=12
   real(pr) :: z0(nc), zi(nc), a
   real(pr) :: z(nc), P, T

   real(pr) :: w(nc), y(nc), x(nc), lnkx(nc), lnky(nc), lnP, lnT, beta
   real(pr) :: XX(2*nc+3), F(2*nc+3), dF(2*nc+3, 2*nc+3)

   real(pr) :: maxP=10000
   integer :: ns, its, i

   write(*, *) test_title("3ph PT envelope")

   ! eos = asphaltenes_srk(z0, zi)
   eos = multicomponent_PR(z0, zi)




   a = 0.7
   z = a * zi + (1-a)*z0


   env = find_hpl(eos, z, t0=500._pr, p0=maxP, max_points=1000)
   write(1, *) env


   P = 0.01
   sat = saturation_temperature(eos, z, P, kind="dew", t0=500._pr)
   env = pt_envelope_2ph(eos, z, sat, points=1000, maximum_pressure=maxP)
   cp = critical_point(eos, z0, zi=0*z0, spec=spec_CP%a, S=0._pr, max_iters=100)
   write(1, *) env

   sat = saturation_pressure(eos, z, T=250._pr, kind="bubble", p0=0.01_pr)
   fl = flash(eos, z=z, T=sat%T, P_spec=sat%P, iters=its)
   env = pt_envelope_2ph(eos, z, sat, points=1000)
   write(1, *) env

   print *, sat%iters
   print *, sat%y

   ! The incipient phase is the unstable gas phase
   w = sat%y
   x = sat%x
   y = fl%y
   beta = 0.5
   
   lnP = log(sat%P)
   lnT = log(sat%T)
   lnkx = log(x/w)
   lnky = log(y/w)
   beta = fl%beta

   ! bub_init: block
   !    !! Correction for fluids that contain asphaltenes, asuming that the
   !    !! asphaltenes are the last component of the composition vector.
   !    real(pr) :: phi_w(nc), phi_y(nc)
   !    y = 0
   !    y(nc) = 1
   !    call eos%lnphi_pt(y, sat%P, sat%T, root_type="liquid", lnphi=phi_y)
   !    call eos%lnphi_pt(w, sat%P, sat%T, root_type="vapor", lnphi=phi_w)
   !    lnKy = phi_w - phi_y
   !    beta = z(nc)
   ! end block bub_init

   print *, "flash"
   print *, fl%iters
   print *, fl%x
   print *, fl%y
   print *, fl%beta
   print *, "flash"

   XX = [lnKx, lnKy, lnP, lnT, beta]

   ns = 2*nc + 2

   print *, "aaaa"
   print "((A,x), *(E13.4,x))", "w0", w
   print "((A,x), *(E13.4,x))", "x0", x
   print "((A,x), *(E13.4,x))", "y0", y
   print *, "aaaa"

   call solve_point(eos, z, ns, XX(ns), XX, F, dF, its, 1000)

   w = z/(XX(2*nc+3)*exp(XX(nc+1:2*nc)) + (1 - XX(2*nc+3))*exp(XX(:nc)))
   x = w * exp(XX(:nc))
   y = w * exp(XX(nc+1:2*nc))
   beta = XX(2*nc+3)

   print *, "w:", w
   print *, "x:", x
   print *, "y:", y
   print *, "beta:", beta
   print * , "P:", exp(XX(2*nc+1))
   print * , "T:", exp(XX(2*nc+2))

   env3 = pt_envelope_3ph(&
      eos, z=z, x0=x, y0=y, w0=w, beta0=beta, &
      P0=exp(XX(2*nc+1)), T0=exp(XX(2*nc+2)), ns0=ns, dS0=0.1_pr, points=10000)

  write(2, *) "ns S T P beta &
               x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 &
               y1 y2 y3 y4 y5 y6 y7 y8 y9 y10 y11 y12 y13 y14 y15 y16 &
               w1 w2 w3 w4 w5 w6 w7 w8 w9 w10 w11 w12 w13 w14 w15 w16"
  do i=1,size(env3%beta)
     write(2, *) env3%ns(i), env3%S(i), env3%T(i), env3%P(i), env3%beta(i), env3%x(i, :), env3%y(i, :), env3%w(i, :)
  end do

end program main
