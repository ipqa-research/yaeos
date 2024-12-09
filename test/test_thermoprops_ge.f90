program main
   use yaeos, only: pr, R
   use yaeos, only: setup_unifac, UNIFAC, Groups
   use auxiliar_functions, only: allclose, rel_error

   implicit none

   integer, parameter :: nc = 3

   type(UNIFAC) :: model
   type(Groups) :: molecules(nc)

   real(pr) :: n(nc), n_aux1(nc), n_aux2(nc), T

   real(pr) :: lngamma(nc), dlngammadT(nc), dlngammadn(nc, nc)
   real(pr) :: lngammadn_num(nc, nc), lngammadT_num(nc), lngamma_aux1(nc)
   real(pr) :: lngamma_aux2(nc), dlngammadn_aux(nc, nc), dlngammadT_aux(nc)

   real(pr) :: Ge

   real(pr) :: He, HeT, Hen(nc)
   real(pr) :: Hen_num(nc), HeT_num, He_aux1, He_aux2, HeT_aux, Hen_aux(nc)

   real(pr) :: Se, SeT, Sen(nc)
   real(pr) :: Sen_num(nc), SeT_num, Se_aux1, Se_aux2, SeT_aux, Sen_aux(nc)

   real(pr) :: dn=0.0001_pr, dT=0.0001_pr

   integer :: i, j

   T = 303.15_pr
   n = [21.48_pr, 11.42_pr, 16.38_pr]

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

   ! Properties on point
   call model%excess_gibbs(n, T, Ge=Ge)

   call model%ln_activity_coefficient(&
      n, T, lngamma=lngamma, dlngammadT=dlngammadT, dlngammadn=dlngammadn &
      )

   call model%excess_enthalpy(n, T, He=He, HeT=HeT, Hen=Hen)

   call model%excess_entropy(n, T, Se=Se, SeT=SeT, Sen=Sen)

   ! ==========================================================================
   ! Obvious test
   ! --------------------------------------------------------------------------
   ! Activity coefficients
   if (.not. rel_error(R * T * sum(n * lngamma), Ge) < 1e-10_pr) then
      error stop "Test Thermoprops GE: ln(gamma) and GE failed"
   end if

   ! Entropy
   if (.not. rel_error(Se, (He - Ge) / T) < 1e-10_pr) then
      error stop "Test Thermoprops GE: Se and GE failed"
   end if

   ! ==========================================================================
   ! Individual calls
   ! --------------------------------------------------------------------------
   ! ln(gamma)
   call model%ln_activity_coefficient(n, T, lngamma=lngamma_aux1)
   call model%ln_activity_coefficient(n, T, dlngammadT=dlngammadT_aux)
   call model%ln_activity_coefficient(n, T, dlngammadn=dlngammadn_aux)

   if (.not. allclose(lngamma, lngamma_aux1, 1.0E-10_pr)) then
      error stop "Test Thermoprops GE: ln(gamma) failed"
   end if

   if (.not. allclose(dlngammadT, dlngammadT_aux, 1.0E-10_pr)) then
      error stop "Test Thermoprops GE: dln(gamma)/dT failed"
   end if

   do i = 1,nc
      if (.not. allclose(dlngammadn(i,:), dlngammadn_aux(i,:),1.0E-10_pr)) then
         error stop "Test Thermoprops GE: dln(gamma)/dn failed"
      end if
   end do

   ! Excess enthalpy
   call model%excess_enthalpy(n, T, He=He_aux1)
   call model%excess_enthalpy(n, T, HeT=HeT_aux)
   call model%excess_enthalpy(n, T, Hen=Hen_aux)

   if (.not. rel_error(He, He_aux1) < 1e-10_pr) then
      error stop "Test Thermoprops GE: He failed"
   end if

   if (.not. rel_error(HeT, HeT_aux) < 1e-10_pr) then
      error stop "Test Thermoprops GE: HeT failed"
   end if

   if (.not. allclose(Hen, Hen_aux, 1.0E-10_pr)) then
      error stop "Test Thermoprops GE: Hen failed"
   end if

   ! Excess entropy
   call model%excess_entropy(n, T, Se=Se_aux1)
   call model%excess_entropy(n, T, SeT=SeT_aux)
   call model%excess_entropy(n, T, Sen=Sen_aux)

   if (.not. rel_error(Se, Se_aux1) < 1e-10_pr) then
      error stop "Test Thermoprops GE: Se failed"
   end if

   if (.not. rel_error(SeT, SeT_aux) < 1e-10_pr) then
      error stop "Test Thermoprops GE: SeT failed"
   end if

   if (.not. allclose(Sen, Sen_aux, 1.0E-10_pr)) then
      error stop "Test Thermoprops GE: Sen failed"
   end if

   ! ==========================================================================
   ! Derivatives
   ! --------------------------------------------------------------------------
   ! dln(gamma)/dT
   call model%ln_activity_coefficient(n, T + dT, lngamma=lngamma_aux1)
   call model%ln_activity_coefficient(n, T - dT, lngamma=lngamma_aux2)

   lngammadT_num = (lngamma_aux1 - lngamma_aux2)/(2.0_pr*dT)

   if (.not. allclose(dlngammadT, lngammadT_num, 1.0E-5_pr)) then
      error stop "Test Thermoprops GE: dln(gamma)/dT failed"
   end if

   ! dln(gamma)/dn
   lngammadn_num = 0.0_pr

   do i=1, nc
      do j=1,nc
         n_aux1 = n
         n_aux2 = n

         n_aux1(j) = n_aux1(j) + dn
         n_aux2(j) = n_aux2(j) - dn

         call model%ln_activity_coefficient(n_aux1, T, lngamma=lngamma_aux1)
         call model%ln_activity_coefficient(n_aux2, T, lngamma=lngamma_aux2)

         lngammadn_num(i,j) = (lngamma_aux1(i) - lngamma_aux2(i))/(2.0_pr * dn)
      end do
   end do

   do i = 1,nc
      if (.not. allclose(dlngammadn(i,:), lngammadn_num(i,:), 1.0E-5_pr)) then
         error stop "Test Thermoprops GE: dln(gamma)/dn failed"
      end if
   end do

   ! dHe/dT
   call model%excess_enthalpy(n, T + dT, He=He_aux1)
   call model%excess_enthalpy(n, T - dT, He=He_aux2)

   HeT_num = (He_aux1 - He_aux2) / (2.0_pr * dT)

   if (.not. rel_error(HeT, HeT_num) < 1e-6) then
      error stop "Test Thermoprops GE: dHe/dT failed"
   end if

   ! dHe/dn
   Hen_num = 0.0_pr

   do i=1, nc
      n_aux1 = n
      n_aux2 = n

      n_aux1(i) = n_aux1(i) + dn
      n_aux2(i) = n_aux2(i) - dn

      call model%excess_enthalpy(n_aux1, T, He=He_aux1)
      call model%excess_enthalpy(n_aux2, T, He=He_aux2)

      Hen_num(i) = (He_aux1 - He_aux2) / (2.0_pr * dn)
   end do

   if (.not. allclose(Hen, Hen_num, 1.0E-6_pr)) then
      error stop "Test Thermoprops GE: dHe/dn failed"
   end if

   ! dSe/dT
   call model%excess_entropy(n, T + dT, Se=Se_aux1)
   call model%excess_entropy(n, T - dT, Se=Se_aux2)

   SeT_num = (Se_aux1 - Se_aux2) / (2.0_pr * dT)

   if (.not. rel_error(SeT, SeT_num) < 1e-6) then
      error stop "Test Thermoprops GE: dSe/dT failed"
   end if

   ! dSe/dn
   Sen_num = 0.0_pr

   do i=1, nc
      n_aux1 = n
      n_aux2 = n

      n_aux1(i) = n_aux1(i) + dn
      n_aux2(i) = n_aux2(i) - dn

      call model%excess_entropy(n_aux1, T, Se=Se_aux1)
      call model%excess_entropy(n_aux2, T, Se=Se_aux2)

      Sen_num(i) = (Se_aux1 - Se_aux2) / (2.0_pr * dn)
   end do

   if (.not. allclose(Sen, Sen_num, 1.0E-6_pr)) then
      error stop "Test Thermoprops GE: dSe/dn failed"
   end if
end program main
