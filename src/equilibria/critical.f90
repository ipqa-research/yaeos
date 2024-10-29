module yaeos__equilibria_critical
   use yaeos__constants, only: pr
   use yaeos__models, only: ArModel
   use yaeos__equilibria_equilibrium_state, only: EquilibriumState

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
      use yaeos__equilibria_equilibrium_state, only: EquilibriumState
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: a0
      real(pr), intent(in) :: z0(:)
      real(pr), intent(in) :: zi(:)
      real(pr), intent(in) :: dS0

      real(pr) :: u(size(z0)), u_new(size(z0))

      real(pr), allocatable :: XS(:, :)
      real(pr) :: X0(3), T, P, V, z(size(z0))

      type(EquilibriumState) :: first

      integer :: i, j, ns, last_point
      ! u = z0
      ! u = zi
      u = (z0 + zi)/sum(z0 + zi)
      z = a0*zi + (1-a0)*z0

      T = sum(model%components%Tc * z)
      P = sum(model%components%Pc * z)
      call model%volume(n=z, P=P, T=T, V=V, root_type="stable")

      X0 = [a0, v, T]
      X0(2:3) = log(X0(2:3))
      ns = 1

      ! call solve_cp(model, X0, ns, X0(ns), z0, zi, u)

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

         ns = maxloc(abs(dXdS), dim=1)
         dS = dXdS(ns)*dS
         dXdS = dXdS/dXdS(ns)

         if (exp(X(2)) < 0.1) then
            ! If the volume is too small, reduce the step size
            do while(abs(dXdS(2)*dS) > abs(0.01 * X(2)))
               dS = dS/2
            end do
         end if
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

   type(EquilibriumState) function critical_point(&
      model, z0, zi, spec, &
      max_iters, u0, V0, T0, a0 &
      )
      use yaeos__math, only: solve_system
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: z0(:)
      real(pr), intent(in) :: zi(:)
      character(len=*), intent(in) :: spec
      integer, intent(in) :: max_iters

      real(pr), optional, intent(in) :: V0
      real(pr), optional, intent(in) :: T0
      real(pr), optional, intent(in) :: a0

      real(pr), optional, intent(in) :: u0(:)

      real(pr) :: X(3)
      integer :: ns
      real(pr) :: S
      real(pr) :: F(3), df(3, 3), dX(3), u(size(z0))
      real(pr) :: V, T, P

      real(pr) :: z(size(z0)), u_new(size(z0)), l
      integer :: i

      ! ========================================================================
      ! Handle the input
      ! ------------------------------------------------------------------------
      if (present(a0)) then
         X(1) = a0
      else
         X(1) = 0.0_pr
      end if

      z = X(1)*zi + (1-X(1))*z0

      if (present(u0)) then
         u = u0
      else
         u = (z0 + zi)/sum(z0 + zi)
      end if

      if (present(T0)) then
         X(3) = log(T0)
      else
         X(3) = log(sum(model%components%Tc * z))
      end if

      if (present(V0)) then
         X(2) = log(V0)
      else
         call model%volume(&
            n=z, P=sum(model%components%Pc * z), T=exp(X(3)), V=X(2), root_type="stable")
         X(2) = log(X(2))
      end if

      select case (spec)
       case("z")
         ns = 1
         S = X(1)
       case("V")
         ns = 2
         S = X(2)
       case("T")
         ns = 3
         S = X(3)
       case default
         stop "Invalid specification"
      end select

      ! ========================================================================
      ! Solve the system of equations
      ! ------------------------------------------------------------------------
      do i=1,max_iters
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

         critical_point%iters = i
      end do

      critical_point%x = z
      critical_point%y = z
      critical_point%Vx = exp(X(2))
      critical_point%Vy = exp(X(2))
      critical_point%T  = exp(X(3))
      call model%pressure(n=z, V=critical_point%Vx, T=critical_point%T, P=critical_point%P)
      critical_point%kind = "critical"

   end function critical_point
end module yaeos__equilibria_critical
