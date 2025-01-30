program main
   use yaeos

   ! Set the variable `model` as a generic `ArModel`
   class(ArModel), allocatable :: model

   ! Set the variables that we're going to use as variable length arrays
   real(pr), allocatable :: n(:), tc(:), pc(:), w(:)

   n = [0.3, 0.7]    ! Number of moles of each component [mol]
   tc = [190, 310]   ! Critical temperatures [K]
   pc = [14, 30]     ! Critical pressures [bar]
   w = [0.001, 0.03] ! Acentric factors [-]

   ! Now we set our model as the PengRobinson76 equation of state.
   model = PengRobinson76(tc, pc, w)

   pressure: block
      real(pr) :: T, T2, dT, V, P, P2, dPdV, dPdT, dPdn(2)

      T = 300.0_pr   ! Set temperature to 300 K
      T2 = 300.1_pr
      dt = T2 - T  
      V = 0.1_pr     ! Set volume to 0.1 L

      ! Calculate pressure and its derivatives
      call model%pressure(n, V, T, P=P, dPdV=dPdV, dPdT=dPdT, dPdn=dPdn)
      call model%pressure(n, V, T2, P=P2)

      print *, "Pressure: ", P
      print *, "dPdV: ", dPdV
      print *, "dPdT: ", dPdT
      print *, "dPdn: ", dPdn

      print *, "numericam: ", (P2 - P) / dt

   end block pressure

end program
