program main
   use yaeos
   type(CriticalLine) :: cl
   type(CubicEoS) :: model
   real(pr) :: zi(2), z0(2), tc(2), pc(2), w(2), eps
   integer :: i

   zi = [0.1_pr, 0.9_pr]

   tc = [304.21, 727.0]
   pc = [73.83, 25.6]
   w =  [0.223621, 0.427556]

   model = PengRobinson76(tc, pc, w)

   eps = 1e-3
   zi = [1-epsilon(eps), eps]
   z0 = [eps, 1 -eps]

   cl = critical_line(model, a0=eps, z0=z0, zi=zi, ds0=1e-3_pr)

   ! print *, size(cl%a)

   do i=1,size(cl%a)
      print *, cl%a(i), cl%T(i), cl%P(i), cl%V(i)
   end do

   ! yaeos.yaeos_c.critical_line(model.id, a0=eps, da0=1e-3, z0=z0, zi=zi, max_points=100)

end program