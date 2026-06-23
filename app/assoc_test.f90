program main
   use yaeos

   implicit none

   integer, parameter :: nc=2
   type(PcSaft) :: pcsaft_adiff

   real :: t0, t1

   real(pr) :: m(nc), sigma(nc), eps_k(nc), kij(nc,nc), n(nc), T, V, z(nc)
   real(pr) :: eps_ab(nc), kap_ab(nc)
   integer :: na(nc), nb(nc)
   integer :: its

   real(pr) :: adAr, adArT, adArV, adArn(nc), adArT2, adArV2, adArTn(nc), adArVn(nc), adArTV, adArn2(nc,nc)
   real(pr) :: alnphi(nc), P


   n = [5.0_pr, 4.0_pr]
   T = 303.15_pr
   V = 100.0_pr

   z = n / sum(n)

   m = [1.3403_pr, 1.0656_pr]
   sigma = [3.8582_pr, 3.0007_pr]
   eps_k = [211.59_pr, 366.51_pr]
   na = [1, 1]
   nb = [1, 1]
   kap_ab = [0.07555_pr, 0.034868_pr]
   eps_ab = [3044.4_pr, 2500.7_pr]

   pcsaft_adiff = init_pcsaft(m, sigma, eps_k, eps_ab=eps_ab, kap_ab=kap_ab, na_sites=na, nb_sites=nb)

   ! call cpu_time(t0)
   ! call pcsaft_adiff%residual_helmholtz(n, V, T, adAr, adArV, adArT, adArTV, adArV2, adArT2, adArn, adArVn, adArTn, adArn2)
   ! call cpu_time(t1)

   call cpu_time(t0)
   call pcsaft_adiff%lnphi_vt(n, V, T, lnPhi=alnphi)
   call pcsaft_adiff%pressure(n, V, T, P)
   call cpu_time(t1)

   print *, "Tiempo addiff: ", t1 - t0
   ! print *, "Ar adiff: ", adAr
   ! print *, "Ar adiff V: ", adArV
   ! print *, "Ar adiff T: ", adArT
   ! print *, "Ar adiff TV: ", adArTV
   ! print *, "Ar adiff V2: ", adArV2
   ! print *, "Ar adiff T2: ", adArT2
   ! print *, "Ar adiff n: ", adArn
   ! print *, "Ar adiff Vn: ", adArVn
   ! print *, "Ar adiff Tn: ", adArTn
   ! print *, "Ar adiff n2 fila 1: ", adArn2(1,:)
   ! print *, "Ar adiff n2 fila 2: ", adArn2(2,:)
   ! print *, "Ar adiff n2 fila 3: ", adArn2(3,:)
   ! print *, "Ar adiff n2 fila 4: ", adArn2(4,:)
   print *, "lnPhi adiff: ", alnphi
   print *, "P adiff: ", P * 100
end program main