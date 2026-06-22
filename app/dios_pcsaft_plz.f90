program main
   use yaeos

   implicit none

   integer, parameter :: nc=4
   type(PcSaft) :: pcsaft_adiff
   type(TPCSAFT) :: pcsaft_tapenade

   real :: t0, t1

   real(pr) :: m(nc), sigma(nc), eps_k(nc), kij(nc,nc), n(nc), T, V, z(nc)
   integer :: its

   real(pr) :: adAr, adArT, adArV, adArn(nc), adArT2, adArV2, adArTn(nc), adArVn(nc), adArTV, adArn2(nc,nc)
   real(pr) :: tAr, tArT, tArV, tArn(nc), tArT2, tArV2, tArTn(nc), tArVn(nc), tArTV, tArn2(nc,nc)
   real(pr) :: alnphi(nc), tlnphi(nc)

   type(EquilibriumState) :: fr


   n = [5.0_pr, 4.0_pr, 10.0_pr, 15.0_pr]
   T = 303.15_pr
   V = 100.0_pr

   z = n / sum(n)

   m = [1.0_pr, 1.6069_pr, 2.002_pr, 2.3316_pr]
   sigma = [3.7039_pr, 3.5206_pr, 3.6184_pr, 3.7086_pr]
   eps_k = [150.03_pr, 191.42_pr, 208.11_pr, 222.88_pr]
   
   kij(1,:) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   kij(2,:) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   kij(3,:) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   kij(4,:) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]

   pcsaft_adiff = init_pcsaft(m, sigma, eps_k, kij)
   pcsaft_tapenade = init_tpcsaft(m, sigma, eps_k, kij)

   print *, "------------------------------------------------------------------"
   print *, "PCSAFT adiff"
   print *, "------------------------------------------------------------------"

   call cpu_time(t0)
   call pcsaft_adiff%residual_helmholtz(n, V, T, adAr, adArV, adArT, adArTV, adArV2, adArT2, adArn, adArVn, adArTn, adArn2)
   call cpu_time(t1)

   call pcsaft_adiff%lnphi_vt(n, V, T, lnPhi=alnphi)

   print *, "Tiempo addiff: ", t1 - t0
   print *, "Ar adiff: ", adAr
   print *, "Ar adiff V: ", adArV
   print *, "Ar adiff T: ", adArT
   print *, "Ar adiff TV: ", adArTV
   print *, "Ar adiff V2: ", adArV2
   print *, "Ar adiff T2: ", adArT2
   print *, "Ar adiff n: ", adArn
   print *, "Ar adiff Vn: ", adArVn
   print *, "Ar adiff Tn: ", adArTn
   print *, "Ar adiff n2 fila 1: ", adArn2(1,:)
   print *, "Ar adiff n2 fila 2: ", adArn2(2,:)
   print *, "Ar adiff n2 fila 3: ", adArn2(3,:)
   print *, "Ar adiff n2 fila 4: ", adArn2(4,:)
   print *, "lnPhi adiff: ", alnphi

   print *, "------------------------------------------------------------------"
   print *, "PCSAFT tapenade"
   print *, "------------------------------------------------------------------"
   call cpu_time(t0)
   call pcsaft_tapenade%residual_helmholtz(n, V, T, tAr, tArV, tArT, tArTV, tArV2, tArT2, tArn, tArVn, tArTn, tArn2)
   call cpu_time(t1)

   call pcsaft_tapenade%lnphi_vt(n, V, T, lnPhi=tlnphi)

   print *, "Tiempo tapenade: ", t1 - t0
   print *, "Ar tapenade: ", tAr
   print *, "Ar tapenade V: ", tArV
   print *, "Ar tapenade T: ", tArT
   print *, "Ar tapenade TV: ", tArTV
   print *, "Ar tapenade V2: ", tArV2
   print *, "Ar tapenade T2: ", tArT2
   print *, "Ar tapenade n: ", tArn
   print *, "Ar tapenade Vn: ", tArVn
   print *, "Ar tapenade Tn: ", tArTn
   print *, "Ar tapenade n2 fila 1: ", tArn2(1,:)
   print *, "Ar tapenade n2 fila 2: ", tArn2(2,:)
   print *, "Ar tapenade n2 fila 3: ", tArn2(3,:)
   print *, "Ar tapenade n2 fila 4: ", tArn2(4,:)
   print *, "lnPhi tapenade: ", tlnphi


   print *, "------------------------------------------------------------------"
   print *, "Flash adiff"
   print *, "------------------------------------------------------------------"
   its = 100

   call cpu_time(t0)
   fr = flash(pcsaft_adiff, z, P_spec=10.0_pr, T=200.0_pr, iters=its)
   call cpu_time(t1)

   print *, "Tiempo Flash adiff: ", t1 - t0
   print *, "Beta", fr%beta

   print *, "------------------------------------------------------------------"
   print *, "Flash tapenade"
   print *, "------------------------------------------------------------------"
   call cpu_time(t0)
   fr = flash(pcsaft_tapenade, z, P_spec=5.0_pr, T=300.0_pr, iters=its)
   call cpu_time(t1)

   print *, "Tiempo Flash tapenade: ", t1 - t0
   print *, "Beta", fr%beta
end program main