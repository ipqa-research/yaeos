program phase_envelope
   use yaeos
   use yaeos__example_tools, only: methane_butane_pr76
   implicit none

   class(ArModel), allocatable :: model ! Model
   type(EquilibriumState) :: bubble ! Saturation point
   type(PTEnvel2) :: envelope

   integer, parameter :: nc=2
   real(pr) :: z(nc)

   ! Methane/ Butane mixture
   model = methane_butane_pr76()

   ! Composition
   z = [0.1_pr, 0.9_pr]

   ! Calculate a bubble point at 100K to initialize the rest of the phase
   ! envelope calculation
   bubble = saturation_pressure(model, z, T=150._pr, kind="bubble")

   ! Calculate the whole phase envelope using the information from the converged
   ! dew point
   envelope = pt_envelope_2ph(model, z, bubble)

   ! Write the resulting phase envelope to the screen
   write(*, *) envelope

   ! Write the resulting phase envelope to a file
   open(1, file="phase_envelope.dat")
      write(1, *) envelope
   close(1)

   ! It is also possible to initialize the phase_envelope with a manually
   ! defined point
   manual: block
      real(pr) :: k(nc), T, P
      T = 150
      P = 1
      k = k_wilson(model, T=T, P=P)

      ! Define manually an initial point
      bubble%kind = "bubble"
      bubble%x = z
      bubble%y = k*z
      bubble%beta = 0

      envelope = pt_envelope_2ph(model, z, bubble)
   end block manual
end program phase_envelope