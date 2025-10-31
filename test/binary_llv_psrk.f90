program main
   use yaeos
   use yaeos__equilibria_binaries
   use yaeos__models_ge_group_contribution_unifac, only: Groups
   use testing_aux, only: test_ok, test_title, assert
   implicit none

   integer, parameter :: nc=2
   type(CubicEoS) :: model
   type(EquilibriumState) :: cp
   type(CriticalLine) :: cl
   type(BinaryThreePhase) :: llv

   real(pr) :: a, V, T, P, z0(nc), zi(nc), z(nc), kij(nc,nc)

   real(pr) :: Tc(nc), Pc(nc), w(nc)
   real(pr) :: lnphi(nc), dlnphidn(nc, nc), lambd(50), as(50)

   real(pr) :: X(7), F(7), dF(7,7), S, F1(7), F2(7)
   real(pr) :: aux(7, 7), XdX(7), eps=1e-5, numdf(7, 7)
   integer :: ns

   integer :: i


   type(Groups) :: molecules(2)

   molecules(1)%groups_ids = [1, 2, 3, 14]
   molecules(1)%number_of_groups = [2, 2, 1, 1]

   molecules(2)%groups_ids = [16]
   molecules(2)%number_of_groups = [1]

   write(*, *) test_title("Three-Phase Line from Critical End Point")

   tc = [577.2, 647.13]
   pc = [39.300000000000004, 220.55]
   w = [0.589027, 0.344861]


   model = PSRK(Tc, Pc, w, molecules)

   a = 0.5

   T = 500
   P = 1990
   z0 = [1, 0]
   zi = [0, 1]

   call find_llcl(model, z0, zi, P, a, V, T)
   z = a * zi + (1-a) * z0

   cp%x = z
   cp%T = T
   cp%P = P
   call model%volume(z, P, T, cp%Vx, root_type="liquid")

   cp = critical_point(model, z0, zi, spec=spec_CP%P, S=log(cp%P), &
      t0=T, a0=a, v0=cp%VX, max_iters=200)
   cl = critical_line(model, a0=a, z0=z0, zi=zi, ns0=4, S0=log(cp%P), &
      ds0=-0.1_pr, first_point=cp, stability_analysis=.true., max_points=3000)
   
   X = log([&
      cl%CEP%x(1)+1e-9, &
      cl%cep%x(1)-1e-9, &
      cl%cep%y(1), &
      cl%cep%Vx+1e-9, &
      cl%cep%Vx-1e-9, &
      cl%cep%Vy, &
      cl%cep%T &
      ])

   ns = 0
   S = exp(X(1)) - exp(X(2))
   
   call three_phase_line_F_solve(model, X, ns, S, F, dF)
   print *, maxval(abs(F))
   llv = binary_llv_from_cep(model, cl%cep)
   i = size(llv%T)
   call assert(llv%T(i) < 100, "Three-phase line reaches low temperatures")
end program main
