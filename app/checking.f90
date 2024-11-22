program main
   use yaeos, only: pr, R
   use yaeos, only: Groups, UNIFAC, setup_unifac
   use yaeos, only: UNIQUAC, setup_uniquac

   integer, parameter :: nc = 3, ng = 4

   type(UNIFAC) :: unif
   type(UNIQUAC) :: uniq
   type(Groups) :: molecules(nc)

   real(pr) :: Ge, Gen(nc), GeT, GeT2, GeTn(nc), Gen2(nc, nc)
   real(pr) :: Ge_i, Gen_i(nc), GeT_i, GeT2_i, GeTn_i(nc), Gen2_i(nc, nc)
   real(pr) :: ln_gammas(nc), b(nc, nc), rs(nc), qs(nc)

   real(pr) :: n(nc), T, n_t

   T = 150
   n = [0.0_pr, 70.0_pr, 10.0_pr]
   n_t = sum(n)

   ! ! Ethane [CH3]
   molecules(1)%groups_ids = [1]
   molecules(1)%number_of_groups = [2]

   ! ! Ethanol [CH3, CH2, OH]
   molecules(2)%groups_ids = [1, 2, 14]
   molecules(2)%number_of_groups = [1, 1, 1]

   ! ! Methylamine [H3C-NH2]
   molecules(3)%groups_ids = [28]
   molecules(3)%number_of_groups = [1]

   unif = setup_unifac(molecules)

   ! ===========================================================================
   ! UNIFAC
   ! ---------------------------------------------------------------------------
   call unif%excess_gibbs(n, T, Ge, GeT, GeT2, Gen, GeTn, Gen2)

   print *, "UNIFAC:"
   print *, "Ge: ", Ge
   print *, "GeT: ", GeT
   print *, "GeT2: ", GeT2
   print *, "Gen: ", Gen
   print *, "GeTn: ", GeTn
   print *, "Gen2: ", Gen2

   ! ===========================================================================
   ! UNIQUAC
   ! ---------------------------------------------------------------------------
   rs = [0.92_pr, 2.1055_pr, 3.1878_pr]
   qs = [1.4_pr, 1.972_pr, 2.4_pr]

   T = 298.15_pr

   ! Calculate bij from DUij. We need -DU/R to get bij
   b(1,:) = [0.0_pr, -526.02_pr, -309.64_pr]
   b(2,:) = [318.06_pr, 0.0_pr, 91.532_pr]
   b(3,:) = [-1325.1_pr, -302.57_pr, 0.0_pr]

   uniq = setup_uniquac(qs, rs, bij=b)

   call uniq%excess_gibbs(n, T, Ge, GeT, GeT2, Gen, GeTn, Gen2)

   print *, "UNIQUAC:"
   print *, "Ge: ", Ge
   print *, "GeT: ", GeT
   print *, "GeT2: ", GeT2
   print *, "Gen: ", Gen
   print *, "GeTn: ", GeTn
   print *, "Gen2: ", Gen2

end program main