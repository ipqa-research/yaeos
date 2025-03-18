program main
   use yaeos
   use fixtures_models, only: multicomponent_PR
   use testing_aux, only: assert, test_title
   implicit none
   
   type(CubicEoS) :: eos
   type(PTEnvel2) :: env
   type(EquilibriumState) :: sat, cp
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
   cp = critical_point(eos, z0, zi=0*z0, spec=spec_CP%a, S=0._pr, max_iters=100)
   call assert(abs(env%cps(1)%T - Tc)/Tc < 1e-1, "Critical Temperature")
   call assert(abs(env%cps(1)%P - Pc)/pc < 1e-1, "Critical Pressure")
end program