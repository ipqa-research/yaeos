module yaeos__equilibria_critical
   use yaeos__constants, only: pr
   use yaeos__models, only: ArModel

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

   type(CriticalLine) function critical_line(model, a0, z0, zi, dS0)
      use yaeos__math_continuation, only: continuation
      use yaeos__math, only: solve_system
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: a0
      real(pr), intent(in) :: z0(:)
      real(pr), intent(in) :: zi(:)
      real(pr), intent(in) :: dS0

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

      call solve_cp(model, X0, ns, X0(ns), z0, zi, u)

      u = [(1, i=1, size(z0))]
      u = u/sum(u)
      u = zi

      XS = continuation(&
         f=foo, X0=X0, ns0=ns, S0=X0(ns), &
         dS0=dS0, max_points=2500, solver_tol=1e-5_pr, &
         update_specification=update_specification &
         )

      last_point = 0
      do i=1, size(XS, 1)
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

      subroutine update_specification(X, ns, S, dS, dXdS, iterations)
         real(pr), intent(in out) :: X(:) !! Vector of variables \(X\)
         integer, intent(in out) :: ns !! Position of specified variable
         real(pr), intent(in out) :: S !! Specification variable value
         real(pr), intent(in out) :: dS !! Step of specification in the method
         real(pr), intent(in out) ::  dXdS(:) !! \(\frac{dX}{dS}\)
         integer, intent(in) :: iterations !! Iterations needed to converge point

         integer :: other(2) = [2,3]

         ns = maxloc(abs(dXdS), dim=1)
         dS = dXdS(ns)*dS
         dXdS = dXdS/dXdS(ns)
      end subroutine update_specification

   end function critical_line

   real(pr) function lambda1(model, X, s, z0, zi, u, u_new)
      use stdlib_linalg, only: eigh, linalg_state_type
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

      type(linalg_state_type) :: stat

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
      
      call eigh(A=M, lambda=lambda, vectors=vectors, err=stat)
      if (.not. stat%ok()) write(*, *) stat%print_msg()
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

      real(pr), parameter :: eps=1e-10_pr

      if(any(zi * X(1) + z0 * (1-X(1))<0) ) return


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