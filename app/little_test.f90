program main
   use yaeos

   implicit none

   real(pr) :: tc(3), pc(3), w(3)
   real(pr) :: n(3), T, P, V, Z
   class(ArModel), allocatable :: model

   real(pr) :: Ar_v, Ar_p
   real(pr) :: ArP, ArT, Arn(3)
   real(pr) :: dT, dP, dn, Arn_n(3)
   real(pr) :: Ar_mas_dx, Ar_menos_dx
   integer :: i

   tc = [190_pr, 310_pr, 400_pr]   ! Critical temperatures [K]
   pc = [14_pr, 30_pr, 50_pr]      ! Critical pressures [bar]
   w = [0.001_pr, 0.03_pr, 0.1_pr] ! Acentric factors [-]

   model = PengRobinson76(tc, pc, w)

   T = 303.15_pr
   P = 1.01325_pr
   n = [1.0_pr, 2.0_pr, 4.0_pr]

   call model%volume(n, P, T, V, root_type="stable")

   z = P*V/(sum(n)*R*T)

   call model%residual_helmholtz(n, V, T, Ar=Ar_v)
   call model%helmholtz_residual_pt(n, P, T, root_type="stable", Ar=Ar_p, ArP=ArP, ArT=ArT, Arn=Arn)

   print *, "Helhol"
   print *, Ar_v - sum(n) * R * T * log(Z)
   print *, Ar_p

   print *, "ArP y ArP numerica"
   dp = 0.00001_pr
   call model%helmholtz_residual_pt(n, P + dp, T, root_type="stable", Ar=Ar_mas_dx)
   call model%helmholtz_residual_pt(n, P - dp, T, root_type="stable", Ar=Ar_menos_dx)

   print *, ArP, (Ar_mas_dx - Ar_menos_dx)/(2.0_pr*dp)

   print *, "ArT y ArT numerica"
   dt = 0.001_pr
   call model%helmholtz_residual_pt(n, P, T + dt, root_type="stable", Ar=Ar_mas_dx)
   call model%helmholtz_residual_pt(n, P, T - dt, root_type="stable", Ar=Ar_menos_dx)

   print *, ArT, (Ar_mas_dx - Ar_menos_dx)/(2.0_pr*dt)

   print *, "Arn y Arn numerica"
   dn = 0.01_pr

   do i = 1, 3
      n(i) = n(i) + dn
      call model%helmholtz_residual_pt(n, P, T, root_type="stable", Ar=Ar_mas_dx)
      n(i) = n(i) - 2.0_pr*dn
      call model%helmholtz_residual_pt(n, P, T, root_type="stable", Ar=Ar_menos_dx)
      n(i) = n(i) + dn
      Arn_n(i) = (Ar_mas_dx - Ar_menos_dx)/(2.0_pr*dn)
   end do

   print *, Arn
   print *, Arn_n
end program main
