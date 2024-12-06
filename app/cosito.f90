program main
   use yaeos, only: pr, R
   use yaeos, only: Groups, setup_unifac, UNIFAC

   type(UNIFAC) :: model

   integer, parameter :: nc = 3, ng = 4

   type(Groups) :: molecules(nc)

   real(pr) :: n(nc), T, n_t

   real(pr) :: lngamma(nc), dlngammadT(nc), dlngammadn(nc, nc)
   real(pr) :: He, HeT, Hen(nc)
   real(pr) :: Se, SeT, Sen(nc)

   T = 150
   n = [20.0, 70.0, 10.0]
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

   ! setup UNIFAC model
   model = setup_unifac(molecules)

   call model%ln_activity_coefficient(n, T, lngamma=lngamma, dlngammadT=dlngammadT, dlngammadn=dlngammadn)
   call model%excess_enthalpy(n, T, He=He, HeT=HeT, Hen=Hen)
   call model%excess_entropy(n, T, Se=Se, SeT=SeT, Sen=Sen)

   print *, 'Activity coefficients:', lngamma, [0.0_pr, 0.0_pr, 0.0_pr]
   print *, 'dln(gamma)/dT:', dlngammadT, [0.0_pr, 0.0_pr, 0.0_pr]
   print *, 'dln(gamma)/dn:', dlngammadn
   
   print *, " "

   print *, 'Excess enthalpy:', He / n_t, -8.12666342863807_pr
   print *, 'Excess entropy:', Se / n_t, -0.03268447167877293_pr

   print *, 'dHe/dT:', HeT / n_t, 0.05391608033744438_pr
   print *, 'dSe/dT:', SeT / n_t, 0.0003594405355829625_pr

   print *, 'dHe/dn:', Hen / n_t, [150.72393613,  -573.71657621, -4412.09526745]
   print *, 'dSe/dn:', Sen / n_t, [-6.01538894_pr, -2.23972167_pr, -4.97564211_pr]


end program main