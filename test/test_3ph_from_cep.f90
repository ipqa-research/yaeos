program main
   use yaeos
   use yaeos__equilibria_binaries
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

   write(*, *) test_title("Three-Phase Line from Critical End Point")

   Tc = [304.21, 540.2]
   Pc = [73.83, 27.4]
   w = [0.223621, 0.349469]

   kij = 0
   kij(1,2) = 0.1
   kij(2,1) = 0.1

   model = PengRobinson76(Tc, Pc, w)

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
   do i=1,size(cl%a)
      write(2, *) cl%a(i), cl%T(i), cl%P(i)
   end do

   print*, cl%CEP

   block
      real(pr) :: Vx, Vy, x1(2)
      x1 = cl%CEP%x
      x1(1) = x1(1) + 1e-7_pr
      x1(2) = 1 - x1(1)
      call model%volume(x1, cl%CEP%P, cl%CEP%T, Vx, root_type="liquid")

      x1 = cl%CEP%x
      x1(1) = x1(1) - 1e-7_pr
      x1(2) = 1 - x1(1)
      call model%volume(x1, cl%CEP%P, cl%CEP%T, Vy, root_type="liquid")

      X = log([&
         cl%CEP%x(1)+1e-7, &
         cl%cep%x(1)-1e-7, &
         cl%cep%y(1), &
         Vx, &
         Vy, &
         cl%cep%Vy, &
         cl%cep%T &
         ])

   end block

   ns = -1
   S = X(4) - X(5)
   ! ns = 0
   ! S = exp(X(2)) - exp(X(1))

   call three_phase_line_F_solve(model, X, ns, S, F, dF)
   print *, maxval(abs(F))
   ! llv = binary_llv_from_cep(model, cl%cep)
   ! i = size(llv%T)

   ! do i=1,size(llv%T)
   !    print *, llv%T(i), llv%P(i)
   ! end do
   ! call assert(llv%T(i) < 100, "Three-phase line reaches low temperatures")


   ! Make numdiff
   ! eps = 1e-3

   ! do i=1,7
   !    XdX = X
   !    XdX(i) = X(i) + eps * X(i)
   !    call three_phase_line_F(model, XdX, ns, S, F1, aux)
   !    XdX(i) = X(i) - eps * X(i)
   !    call three_phase_line_F(model, XdX, ns, S, F2, aux)
   !    numdF(:,i) = (F1 - F2) / (2 * eps * X(i))

   !    print *, i
   !    print "(*(E14.4,x))", dF(:, i)
   !    print "(*(E14.4,x))", dF(:, i)
   ! end do
   ! print *, maxval(abs(dF - numdF)/abs(numdF)), maxloc(abs(dF - numdF))
   ! call three_phase_line_F_solve(model, X, ns, S, F, dF)
   ! print *, exp(X)
   ! print *, maxval(abs(F))





end program main
