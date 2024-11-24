program main
   use yaeos
   use fixtures_models, only: binary_PR76
   implicit none

   type(CubicEoS) :: model
   type(PTEnvel2) :: env
   real(pr) :: z(2), T, P

   model = binary_PR76()

   z = [0.5_pr, 0.5_pr]

   env = find_hpl(model, z, T0=500._pr, p0=150._pr)

   if(abs(env%points(1)%T - 152.55) > 0.01) error stop "HPLL Failed T1"
   if(abs(env%points(1)%P - 150) > 0.1) error stop "HPLL Failed T1"
   if (abs(env%points(size(env%points))%P) > 0.1) error stop "HPLL Failed Pend"

end program