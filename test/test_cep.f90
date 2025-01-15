program main
   use yaeos
   use yaeos__equilibria_critical, only: f_cep, f_cep2, df_cep
   use fixtures_models, only: binary_PR76
   implicit none

   integer, parameter :: nc=2

   type(CubicEoS) :: model
   type(EquilibriumState) :: cp
   real(pr) :: X(nc+4)
   real(pr) :: a = 0.3
   real(pr) :: z0(nc) = [1, 0]
   real(pr) :: zi(nc) = [0, 1]
   real(pr) :: z(nc)

   integer :: ns=2
   real(pr) :: s = 0.5
   real(pr) :: u(nc) = 0.5

   real(pr) :: f(nc+4), df(nc+4, nc+4), df_num(nc+4, nc+4), dx=1e-7
   integer :: i
   model = binary_PR76()
   
   cp = critical_point(model, z0, zi, spec=spec_CP%a, S=a, max_iters=300)

   X = [0.5_pr, 0.5_pr, log(cp%Vx), log(0.1_pr), log(250._pr), 0.5_pr]

   f = f_cep2(model, X, ns, s, z0, zi, u)
   df = df_cep(model, X, ns, s, z0, zi, u)

end program