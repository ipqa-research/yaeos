program phase_diagram
   use yaeos, only: &
      ArModel, saturation_pressure, saturation_temperature, EquilibriaState, PengRobinson76, pr, &
      pt_envelope_2ph, PTEnvel2
   implicit none

   class(ArModel), allocatable :: model
   type(EquilibriaState) :: bubble, dew
   type(PTEnvel2) :: envelope

   integer, parameter :: nc=2
   integer :: i
   real(pr) :: n(nc), tc(nc), pc(nc), w(nc), kij(nc, nc), lij(nc, nc), T

   ! Methane/ Butane mixture
   n = [0.1_pr, 0.9_pr]       ! Composition
   tc = [190.564, 425.12]     ! Critical temperatures
   pc = [45.99, 37.96]        ! Critical pressures
   w = [0.0115478, 0.200164]  ! Acentric factors

   kij = 0
   kij(1, 2) = 0.2
   kij(2, 1) = kij(1, 2)

   lij = kij / 2 

   model = PengRobinson76(tc, pc, w, kij, lij)

   ! Calculate a dew point
   T = 100
   bubble = saturation_pressure(model, n, T, kind="bubble", p0=40._pr)
   print *, bubble%T, bubble%p, bubble%iters

   ! Calculate the whole phase envelope using the information from the converged
   ! dew point
   envelope = pt_envelope_2ph(&
      model, n, bubble%y, T0=bubble%T, P0=bubble%p, delta_0=0.01_pr&
   )

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
