program main
use yaeos
   implicit none
   integer, parameter :: nc=15
   integer :: i, j
   integer, allocatable :: idx(:)

   real(pr) :: z(nc), tc(nc), pc(nc), acen(nc)

   real(pr) :: w(nc), mintpd
   real(pr) :: all_minima(nc,  nc+1)

   real(pr) :: P, T
   type(CubicEoS) :: model

   z=[0.0048, 0.00919, 0.43391, 0.1101, 0.06544, 0.00789, 0.03787, 0.01279, 0.02248, 0.02698, 0.22738, 0.03747, 0.0023, 0.00054, 0.00086]
   Tc=[126.2, 304.2, 190.6, 305.4, 369.8, 408.1, 425.2, 460.4, 469.6, 507.4, 691.81, 956.77, 1118.6, 1325.03, 1445.73]
   Pc=[33.94, 73.76, 46., 48.84, 42.46, 36.48, 38., 33.84, 33.74, 29.69, 19.46, 13.08, 10.66, 10.28, 17.3]
   acen=[0.04, 0.225, 0.008, 0.098, 0.152, 0.176, 0.193, 0.227, 0.251, 0.296, 0.68, 1.208, 0.949, 0.182, 1.274]

   model = PengRobinson76(Tc, Pc, acen)


   P = 40
   T = 250
   
   call min_tpd(model, z, P, T, mintpd, w, all_minima)
   
   z = all_minima(15, :nc)
   call min_tpd(model, z, P, T, mintpd, w, all_minima)
   
   print "(A,2x,*(E15.6,x))", "z: ", z
   print "(A,2x,*(E15.6,x))", "Minimal tm: ", mintpd
   print "(A,2x,*(E15.6,x))", "Composition: ", w
   print *, "=================================================================="
   do i=1,nc
      print "(I2,x,L,*(E15.6,x))", i, sum(abs((all_minima(i, :nc) - z)/z)) < 1e-3, all_minima(i, :nc), all_minima(i, nc+1)
   end do
end program