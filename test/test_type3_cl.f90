program main
   use yaeos
   use testing_aux, only: test_title, assert
   implicit none
   integer, parameter :: nc = 2
   type(CubicEoS) :: eos
   real(pr) :: z0(nc), zi(nc), TC(nc), PC(nc), w(nc), kij(nc, nc)
   real(pr) :: T, P, S
   type(CriticalLine) :: cl
   integer :: last


   write(*, *) test_title("DEJA VU")
   
   z0 = [0.0_pr, 1.0_pr]
   zi = [1.0_pr, 0.0_pr]

   Tc = [304.21, 693.0]
   Pc = [73.83, 15.7]
   w =  [0.223621, 0.643017]
   kij = reshape([0.0_pr, 0.08_pr, 0.08_pr, 0.0_pr], [nc, nc])
   eos = PengRobinson76(Tc, Pc, w, kij=kij)

   cl = critical_line(eos, a0=1e-2_pr, z0=z0, zi=zi, ns0=spec_CP%a, S0=1e-2_pr, ds0=1e-1_pr)

   last = size(cl%T)

   call assert(cl%P(last) > 2000, "Critical line pressure should be greater than 2000 bar")


end program main
