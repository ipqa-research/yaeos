!! Calculate the critical line of the binary system of N2-nC7. Starting from
!! N2
program main
   use yaeos
   use testing_aux, only: test_title, assert
   implicit none
   integer, parameter :: nc=2
   type(CubicEoS) :: model
   type(CriticalLine) :: cl
   type(EquilibriumState) :: cp
   real(pr) :: Tc(nc), Pc(nc), w(nc), z0(nc), zi(nc)
   integer :: i

   write(*, *) test_title("CEP N2-nC7")

   Tc = [126.2, 540.2]
   Pc = [34., 27.4]
   w = [0.0377, 0.3495]

   model = PengRobinson76(Tc, Pc, w)
   z0 = [0, 1]
   zi = [1, 0]
   cl = critical_line(&
      model, z0=z0, zi=zi, ns0=1, s0=0.9999_pr, a0=0.9999_pr, ds0=-1e-7_pr, &
      stability_analysis=.true., max_points=100)

   call assert(abs(cl%CEP%T - 126.2) < 0.1, "Critical end point T")
   call assert(abs(cl%CEP%P - 34.) < 0.1, "Critical end point P")
   call assert(abs(cl%CEP%x(1) - 0.9999) < 1e-3, "Critical end point x(1)")
   call assert(abs(cl%CEP%y(1) - 0.345) < 1e-3, "Critical end point y(1)")
end program main
