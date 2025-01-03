program main
   use yaeos
   use fixtures_models, only: multicomponent_PR
   use testing_aux, only: assert, test_title
   implicit none
   
   type(CubicEoS) :: eos
   type(PTEnvel2) :: env
   type(EquilibriumState) :: sat
   integer, parameter :: nc=12
   real(pr) :: z0(nc), zi(nc)
   real(pr) :: z(nc), P, T

   real(pr) :: Tc=699.059
   real(pr) :: Pc=180.226


   write(*, *) test_title("PT envelope test multicomponent")

   eos = multicomponent_PR(z0, zi)

   z = z0
   P = 0.0001
   sat = saturation_temperature(eos, z, P, kind="dew")

   env = pt_envelope_2ph(eos, z, sat)

   call assert(abs(env%cps(1)%T - Tc) < 1e-3, "Critical Temperature")
   call assert(abs(env%cps(1)%P - Pc) < 1e-3, "Critical Pressure")
   write(1, *) env
end program