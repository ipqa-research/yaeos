program phase_diagram
   use yaeos, only: ArModel, saturation_pressure, EquilibriaState, PengRobinson76, pr
   use yaeos__phase_equilibria_boundaries_phase_envelopes_pt, only: pt_boundary_2ph, PTEnvel2
   implicit none

   class(ArModel), allocatable :: model
   type(EquilibriaState) :: bubble, dew
   type(PTEnvel2) :: envelope

   integer, parameter :: nc=2
   integer :: i
   real(pr) :: n(nc), tc(nc), pc(nc), w(nc)

   ! Methane/ Butane mixture
   n = [0.4, 0.6]                      ! Composition
   tc = [190.564, 617.7]              ! Critical temperatures
   pc = [45.99, 21.1]                 ! Critical pressures
   w = [0.0115478, 0.492328]           ! Acentric factors

   ! Use the PengRobinson76 model
   model = PengRobinson76(tc, pc, w)

   ! Calculate a bubble point
   bubble = saturation_pressure(model, n, 150._pr, kind="bubble")

   ! Calculate the whole phase envelope using the bubble point as an initial
   ! state
   envelope = pt_boundary_2ph(model, n, bubble)

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