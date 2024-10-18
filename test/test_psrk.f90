program main
   !! Running a PSRK cubic equation of state example/test
   use yaeos

   implicit none

   ! ===========================================================================
   ! Definition of variables that will be used
   ! ---------------------------------------------------------------------------
   type(CubicEoS) :: eos
   type(Groups) :: molecules(2)
   type(EquilibriumState) :: sat
   type(PTEnvel2) :: env
   real(pr) :: tc(2), pc(2), w(2), n(2)
   real(pr) :: C(2, 3)

   real(pr) :: v, lnphi(2)

   real(pr) :: pressures(7) = [40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]
   real(pr) :: temperatures(7) = [450.0, 460.0, 470.0, 480.0, 490.0, 500.0, 510.0]

   integer :: i

   ! ===========================================================================
   ! Definition of required parameters 
   ! ---------------------------------------------------------------------------

   ! Critical constants
   tc = [304.21_pr, 553.8_pr]
   pc = [7.383e6_pr, 4.080358e6_pr] / 1e5
   w = [0.223621_pr, 0.213_pr]

   ! Molecules groups
   molecules(1)%groups_ids = [117]
   molecules(1)%number_of_groups = [1]
   
   molecules(2)%groups_ids = [2]
   molecules(2)%number_of_groups = [6]

   ! Mathias-Copeman parameters
   C(1, :) = [0.8255_pr, 0.16755_pr, -1.7039_pr]
   C(2, :) = [0.84082_pr, -0.39847_pr, 0.94148_pr]


   ! Setting up the equation of state
   eos = PSRK(&
      tc=tc, &
      pc=pc, &
      w=w, &
      molecules=molecules, &
      c1=C(:, 1), c2=C(:, 2), c3=C(:, 3) &
   )

   ! Molar numbers of each component
   n = [60, 40]

   do i=1,7
      call eos%lnphi_pt(n=n, P=pressures(i), T=temperatures(i), lnphi=lnphi, root_type="stable")
      call eos%volume(n=n, P=pressures(i), T=temperatures(i), V=V, root_type="stable")
      sat = saturation_pressure(eos, n, temperatures(i), kind="bubble", p0=130._pr)

      print *, pressures(i), temperatures(i), exp(lnphi)
      print *, pressures(i), temperatures(i), V
   end do
end program main
