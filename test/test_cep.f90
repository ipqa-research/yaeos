program main
   use yaeos
   use yaeos__equilibria_binaries
   use testing_aux, only: test_ok, test_title, assert
   implicit none

   integer, parameter :: nc=2
   type(CubicEoS) :: model
   type(EquilibriumState) :: cp
   type(CriticalLine) :: cl

   real(pr) :: a, V, T, P, z0(nc), zi(nc), z(nc), aux, kij(nc,nc)

   real(pr) :: Tc(nc), Pc(nc), w(nc)
   real(pr) :: lnphi(nc), dlnphidn(nc, nc), lambd(50), as(50)

   integer :: i, tries


   write(*, *) test_title("CRITICAL END POINTS")

   Tc = [304.21, 540.2]
   Pc = [73.83, 27.4]
   w = [0.223621, 0.349469]

   kij = 0
   kij(1,2) = 0.1
   kij(2,1) = 0.1

   model = PengRobinson76(Tc, Pc, w)

   a = 0.5

   T = 500
   P = 1990
   z0 = [1, 0]
   zi = [0, 1]

   call find_llcl(model, z0, zi, P, a, V, T)
   z = a * zi + (1-a) * z0

   cp%x = z
   cp%T = T
   cp%P = P
   call model%volume(z, P, T, aux, root_type="liquid")
   cp%Vx = aux

   cp = critical_point(model, z0, zi, spec=spec_CP%P, S=log(cp%P), &
      t0=T, a0=a, v0=cp%VX, max_iters=200)
   cl = critical_line(model, a0=a, z0=z0, zi=zi, ns0=4, S0=log(cp%P), &
      ds0=-0.1_pr, first_point=cp, stability_analysis=.true., max_points=3000)

   call assert(abs(cl%CEP%T - 134.79) < 0.1, "Critical end point T")
   call assert(abs(cl%CEP%P - 1.21e-2) < 0.01, "Critical end point P")

end program main
