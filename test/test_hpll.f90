program main
   use yaeos
   use fixtures_models, only: binary_PR76
   implicit none

   type(CubicEoS) :: model
   type(PTEnvel2) :: env
   real(pr) :: z(2), T, P

   integer :: i

   model = binary_PR76()

   z = [0.5_pr, 0.5_pr]

   env = find_hpl(model, z, T0=500._pr, p0=1500._pr, max_points=1000)

   if(abs(env%points(1)%T - 103.69) > 0.01) error stop "HPLL Failed T1"
   if(abs(env%points(1)%P - 1500) > 0.1) error stop "HPLL Failed T1"
   if (abs(env%points(size(env%points))%P) > 0.15) error stop "HPLL Failed Pend"

end program main