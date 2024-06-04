program phase_envelope
   use yaeos, only: pr, ArModel, EquilibriaState, PTEnvel2, PengRobinson76, saturation_pressure, pt_envelope_2ph
   implicit none

   class(ArModel), allocatable :: model ! model
   type(EquilibriaState) :: bubble !
   type(PTEnvel2) :: envelope

   integer, parameter :: nc=2
   integer :: i
   real(pr) :: z(nc)

   ! Methane/ Butane mixture
   model = PengRobinson76( &
      [305.32_pr, 540.2_pr], &
      [48.72_pr, 27.4_pr], &
      [0.099493_pr, 0.349469_pr] &
      )

   ! Composition
   z = [0.1_pr, 0.9_pr]

   ! Calculate a bubble point at 100K to initialize the rest of the phase
   ! envelope calculation
   bubble = saturation_pressure(model, z, T=230._pr, kind="bubble", p0=5.0_pr)

   ! Calculate the whole phase envelope using the information from the converged
   ! dew point
   envelope = pt_envelope_2ph(model, z, bubble)

   ! Write the resulting phase envelope to the screen
   do i = 1, size(envelope%points)
      print *, envelope%points(i)%kind
   end do
end program phase_envelope
