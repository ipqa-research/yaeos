program main
   use yaeos
   use testing_aux, only: test_title, assert
   use fixtures_models, only: binary_PR76
   implicit none

   type(CubicEoS) :: model
   type(PTEnvel2) :: env
   real(pr) :: z(2), T, P

   integer :: i

   write(*, *) test_title("ISOLATED HPLL LINE")

   model = binary_PR76()

   z = [0.5_pr, 0.5_pr]

   env = find_hpl(model, z, T0=500._pr, p0=1500._pr, max_points=1000)

   call assert(abs(env%points(1)%T - 103.69) < 0.01, "HPLL Failed T1")
   call assert(abs(env%points(1)%P - 1500) < 0.1, "HPLL Failed P1")
   call assert(abs(env%points(size(env%points))%P) < 0.15, "HPLL Failed Pend")

end program main