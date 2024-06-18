program main
   use yaeos, only: R, pr, excess_gibbs, Groups, setup_unifac, UNIFAC

   type(UNIFAC) :: model

   integer, parameter :: nc = 3, ng = 4

   type(Groups) :: molecules(nc)

   real(pr) :: Ge, Gen(nc), GeT, GeT2, GeTn(nc), Gen2(nc, nc)

   real(pr) :: n(nc), ln_gammas(nc), T

   T = 150.0_pr
   n = [2.0_pr, 7.0_pr, 1.0_pr]

   ! Ethane [CH3]
   molecules(1)%groups_ids = [1]
   molecules(1)%number_of_groups = [2]

   ! Ethanol [CH3, CH2, OH]
   molecules(2)%groups_ids = [1, 2, 14]
   molecules(2)%number_of_groups = [1, 1, 1]

   ! Methylamine [H3C-NH2]
   molecules(3)%groups_ids = [28]
   molecules(3)%number_of_groups = [1]

   ! setup UNIFAC model
   model = setup_unifac(molecules)

   ! Call all Ge and derivatives
   call excess_gibbs(model, n, T, Ge, GeT, GeT2, Gen, GeTn, Gen2)

   print *, "Ge: ", Ge
   print *, "GeT: ", GeT
   print *, "GeT2: ", GeT2
   print *, "Gen: ", Gen
   print *, "GeTn: ", GeTn
   print *, "Gen2:"
   print *, Gen2(1,:)
   print *, Gen2(2,:)
   print *, Gen2(3,:)

   print *, "ln_gammas: ", Gen / R / T

   call model%ln_activity_coefficient(n, T, ln_gammas)
   print *, "ln_gammas: ", ln_gammas
end program main
