program phase_diagram
   use yaeos, only: &
      ArModel, saturation_pressure, EquilibriaState, PengRobinson76, pr, &
      pt_envelope_2ph, PTEnvel2
   implicit none

   class(ArModel), allocatable :: model
   type(EquilibriaState) :: bubble, dew
   type(PTEnvel2) :: envelope

   integer, parameter :: nc=2
   integer :: i
   real(pr) :: n(nc), tc(nc), pc(nc), w(nc), kij(nc, nc), lij(nc, nc)

   ! Methane/ Butane mixture
   n = [0.4, 0.6]                      ! Composition
   tc = [190.564, 425.12]              ! Critical temperatures
   pc = [45.99, 37.96]                 ! Critical pressures
   w = [0.0115478, 0.200164]           ! Acentric factors

   kij = reshape([0., 0.1, 0.1, 0.], [nc,nc]) 
   lij = kij / 2 
   model = PengRobinson76(tc, pc, w, kij, lij)

   ! Calculate a bubble point
   bubble = saturation_pressure(model, n, 150._pr, kind="bubble")
   print *, "# Initial bubble point"
   print *, bubble%t, bubble%p

   ! Calculate the whole phase envelope using the bubble point as an initial
   ! state
   envelope = pt_envelope_2ph(model, n, bubble)

   print *, "# Envelope"
   print *, "T              P"
   do i=1,size(envelope%points)
      print *, envelope%points(i)%T, envelope%points(i)%P
   end do

   print *, ""
   print *, ""

   print *, "# Critical"
   print *, "T              P"
   do i=1,size(envelope%cps)
      print *, envelope%cps(i)%T,  envelope%cps(i)%P
   end do
end program phase_diagram