program phase_envelope
   use yaeos
   use yaeos__example_tools, only: methane_butane_pr76
   implicit none

   class(ArModel), allocatable :: model ! model
   type(EquilibriaState) :: bubble ! 
   type(PTEnvel2) :: envelope

   integer, parameter :: nc=2
   real(pr) :: z(nc)

   ! Methane/ Butane mixture
   model = methane_butane_pr76()

   ! Composition
   z = [0.1_pr, 0.9_pr]

   ! Calculate a bubble point at 100K to initialize the rest of the phase
   ! envelope calculation
   bubble = saturation_pressure(model, z, T=100._pr, kind="bubble", p0=40._pr)

   ! Calculate the whole phase envelope using the information from the converged
   ! dew point
   envelope = pt_envelope_2ph(model, z, bubble)

   ! Write the resulting phase envelope to the screen
   write(*, *) envelope
end program phase_envelope