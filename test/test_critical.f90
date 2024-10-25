module yaeos__equilibria_critical
   use yaeos__constants, only: pr
   use yaeos__models, only: ArModel
   use stdlib_linalg, only: eigh

   implicit none

   type :: CriticalLine
      real(pr), allocatable :: a(:)
      real(pr), allocatable :: z0(:)
      real(pr), allocatable :: zi(:)
      real(pr), allocatable :: P(:)
      real(pr), allocatable :: V(:)
      real(pr), allocatable :: T(:)
   end type CriticalLine

contains

   type(CriticalLine) function critical_line(model, a0, z0, zi)
      use yaeos__math_continuation, only: continuation
      use yaeos__math, only: solve_system
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: a0
      real(pr), intent(in) :: z0(:)
      real(pr), intent(in) :: zi(:)

      real(pr) :: u(size(z0)), u_new(size(z0))

      real(pr), allocatable :: XS(:, :)
      real(pr) :: X0(3), T, P, V, z(size(z0))

      integer :: i, j, ns, last_point
      u = z0

      T = sum(model%components%Tc * z0)
      P = sum(model%components%Pc * z0)
      call model%volume(n=z0, P=P, T=T, V=V, root_type="stable")

      X0 = [a0, v, T]
      X0(2:3) = log(X0(2:3))
      ns = 1

      u = [(1, i=1, size(z0))]
      u = u/sum(u)

      XS = continuation(&
         f=foo, X0=X0, ns0=ns, S0=X0(ns), &
         dS0=0.1_pr, max_points=50, solver_tol=1e-5_pr &
         )

      last_point = 0
      do i=1, size(XS, 1)
         print *, XS(i, 1), exp(XS(i, 2:3))
         if (all(abs(XS(i, :)) < 0.001)) exit
         last_point = i + 1
      end do

      XS = XS(1:last_point-1, :)

      critical_line%z0 = z0
      critical_line%zi = zi
      critical_line%a =  XS(:, 1)
      critical_line%V = exp(XS(:, 2))
      critical_line%T = exp(XS(:, 3))

      allocate(critical_line%P(size(critical_line%a)))
      do i=1, size(critical_line%a)
         z = critical_line%a(i)*zi + (1-critical_line%a(i))*z0

         call model%pressure(&
            n=z, V=critical_line%V(i), T=critical_line%T(i), &
            P=critical_line%P(i))
      end do

   contains

      subroutine foo(X, ns, S, F, dF, dFdS)
         real(pr), intent(in) :: X(:)
         integer, intent(in) :: ns
         real(pr), intent(in) :: S
         real(pr), intent(out) :: F(:)
         real(pr), intent(out) :: dF(:, :)
         real(pr), intent(out) :: dFdS(:)
         real(pr) :: l1

         real(pr) :: z(size(u)), Xsol(3)

         if (X(1) > 1) then
            return
         end if

         F = F_critical(model, X, ns, S, z0, zi, u)
         df = df_critical(model, X, ns, S, z0, zi, u)
         l1 = lambda1(model=model, X=X, s=0.0_pr, z0=z0, zi=zi, u=u, u_new=u_new)
         u = u_new
         dFdS = 0
         dFdS(size(dFdS)) = -1
      end subroutine foo
   end function critical_line


   real(pr) function lambda1(model, X, s, z0, zi, u, u_new)
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: z0(:)
      real(pr), intent(in) :: zi(:)
      real(pr), intent(in) :: s
      real(pr), intent(in) :: X(3)
      real(pr), intent(in) :: u(:)
      real(pr), optional, intent(out) :: u_new(:)

      real(pr) :: n(size(z0)), V, T
      real(pr) :: dlnf_dn(size(z0), size(z0))
      real(pr) :: lambda(size(z0)), vectors(size(z0), size(z0))

      integer :: i, j, nc
      real(pr) :: M(size(z0), size(z0)), z(size(z0))

      nc = size(z0)

      z = X(1) * zi + (1-X(1))*z0
      n = z + s * u * sqrt(z)
      V = exp(X(2))
      T = exp(X(3))

      call model%lnfug_vt(n=n, V=V, T=T, dlnfdn=dlnf_dn)

      do i=1,nc
         do j=1,nc
            M(i, j) = sqrt(z(i)*z(j)) * dlnf_dn(i, j)
         end do
      end do

      call eigh(A=M, lambda=lambda, vectors=vectors)
      lambda1 = minval(lambda)
      if (present(u_new)) u_new = vectors(:, minloc(lambda, dim=1))
   end function lambda1

   function F_critical(model, X, ns, S, z0, zi, u)
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: X(3)
      integer, intent(in) :: ns
      real(pr), intent(in) :: S
      real(pr), intent(in) :: z0(:)
      real(pr), intent(in) :: zi(:)
      real(pr), intent(in) :: u(:)

      real(pr) :: F_critical(3)
      real(pr) :: z(size(u))

      real(pr), parameter :: eps=1e-5_pr

      F_critical(1) = lambda1(model=model, X=X, s=0.0_pr, z0=z0, zi=zi, u=u)
      F_critical(2) = (&
         lambda1(model=model, X=X, s= eps, zi=zi, z0=z0, u=u) &
         - lambda1(model=model, X=X, s=-eps, zi=zi, z0=z0, u=u))/(2*eps)
      F_critical(3) = X(ns) - S
   end function F_critical

   function df_critical(model, X, ns, S, z0, zi, u)
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: X(3)
      integer, intent(in) :: ns
      real(pr), intent(in) :: S
      real(pr), intent(in) :: z0(:)
      real(pr), intent(in) :: zi(:)
      real(pr), intent(in) :: u(:)
      real(pr) :: df_critical(3, 3)

      real(pr), parameter :: eps=1e-5_pr

      real(pr) :: dx(3), F1(3), F2(3)

      integer :: i

      df_critical = 0
      do i=1,3
         dx = 0
         dx(i) = eps
         F2 = F_critical(model, X+dx, ns, S, z0, zi, u)
         F1 = F_critical(model, X-dx, ns, S, z0, zi, u)
         df_critical(:, i) = (F2 - F1)/(2*eps)
      end do
      ! df_critical(3, :) = 0
      ! df_critical(3, ns) = 1
   end function df_critical

   subroutine solve_cp(model, X, ns, S, z0, zi, u)
      use yaeos__math, only: solve_system
      class(ArModel), intent(in) :: model
      real(pr), intent(inout) :: X(3)
      integer, intent(in) :: ns
      real(pr), intent(in) :: S
      real(pr), intent(in) :: z0(:)
      real(pr), intent(in) :: zi(:)
      real(pr), intent(in out) :: u(:)

      real(pr) :: F(3), df(3, 3), dX(3)

      real(pr) :: z(size(z0)), u_new(size(z0)), l
      integer :: i

      do i=1,250
         F = F_critical(model, X, ns, S, z0, zi, u)
         df = df_critical(model, X, ns, S, z0, zi, u)
         dX = solve_system(A=df, b=-F)

         do while(maxval(abs(dX/X)) > 1e-1)
            dX = dX/10
         end do

         if (maxval(abs(X)) < 1e-5) exit

         X = X + dX
         l = lambda1(model, X, 0.0_pr, z0, zi, u, u_new)
         u = u_new
      end do

   end subroutine solve_cp
end module yaeos__equilibria_critical

program main
   use yaeos
   use yaeos__math, only: solve_system
   use stdlib_linalg, only: eigh
   use fixtures_models, only: binary_PR76
   use yaeos__equilibria_critical, only: &
      lambda1, F_critical, df_critical, CriticalLine, critical_line, solve_cp
   implicit none

   integer, parameter :: nc=2

   type(CubicEoS) :: model
   type(EquilibriumState) :: sat
   type(PTEnvel2) :: env

   type(CriticalLine) :: cl

   real(pr) ::  V, T, P, a

   real(pr) :: z(nc)
   real(pr) :: z0(nc)
   real(pr) :: zi(nc)

   real(pr) :: u(nc)
   integer :: ns
   real(pr) :: S

   real(pr) :: F(3), X(3)
   integer :: i, j

   model = binary_PR76()
   z0 = [0, 1]
   zi = [1, 0]
   u = [1, 0]

   a = real(1, pr)/100._pr

   cl = critical_line(model, a, z0, zi)
   do i=1, size(cl%a)
      write(2, *) cl%a(i), cl%V(i), cl%T(i), cl%P(i)
   end do

   z = a*zi + (1-a)*z0
   T = sum(model%components%Tc * z)
   P = sum(model%components%Pc * z)
   call model%volume(n=z, P=P, T=T, V=V, root_type="stable")
   X = [a, log(V), log(T)]
   ns = 1
   S = X(ns)

   do j=1, 99, 5
      a = real(j, pr)/100
      z = a*zi + (1-a)*z0
      sat = saturation_pressure(model, z, T=150._pr, kind="bubble")
      env = pt_envelope_2ph(model, z, sat)
      X = [a, X(2), X(3)]
      S = X(ns)

      call solve_cp(model, X, ns, S, z0, zi, u)
      F = F_critical(model, x, ns, s, z0, zi, u)
      call model%pressure(n=z, V=exp(X(2)), T=exp(X(3)), P=P)
      write(11, *) X(1), exp(X(2)), exp(X(3)), P
      write(1, *) a, env%cps
   end do

contains

   type(CubicEoS) function get_model()
      use yaeos__models, only: SoaveRedlichKwong
      real(pr) :: tc(3), pc(3), w(3)
      Tc=  [190.564, 304.21, 617.7]
      Pc=  [45.99, 73.83000000000001, 21.1]
      w=  [0.0115478, 0.223621, 0.492328]
      ! Vc=  [0.09859999999999998, 0.09399999999999999, 0.5999999999999999]
      get_model = SoaveRedlichKwong(tc, pc, w)
   end function get_model
end program main
