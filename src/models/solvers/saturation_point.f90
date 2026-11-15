module nonlinear_solvers
   !! # Damped Newton solver
   !! Solve small nonlinear systems using Newton steps and Armijo backtracking.
   !!
   !! # Description
   !! `newton_solve` accepts a callback returning residuals and a dense Jacobian,
   !! updates the initial guess in place, and reports a convergence code. It is
   !! imported by `solve_TP` in this source file, although that routine
   !! currently calls `newton` from `yaeos__math_nonlinearsolvers` instead.
   !!
   !! # Examples
   !! ```fortran
   !! call newton_solve(fun, x, tol, max_iter, its, info)
   !! ```
   use yaeos__constants, only: dp => pr
   use yaeos__math_linalg, only: solve_system
   implicit none
contains

   !============================================================
   !  Damped Newton solver with Armijo backtracking
   !  Single callback: fun(x, f, J)
   !============================================================
   subroutine newton_solve(fun, x, tol, max_iter, its, info)
      !! # newton_solve
      !! Solve a nonlinear system by damped Newton iteration.
      !!
      !! # Description
      !! Forms the Newton direction from \(J\,dx=-F\), then applies
      !! backtracking to reduce the squared-residual merit function. `info=0`
      !! indicates residual convergence (max-norm below `tol`); `info=1`
      !! indicates the iteration limit.
      !!
      !! # Examples
      !! See the module-level callback example.
      implicit none

      ! ------ Interface to user-supplied function ------
      interface
         subroutine fun(x, f, J)
            import dp
            real(dp), intent(in)  :: x(:) !! Current unknown vector.
            real(dp), intent(out) :: f(:) !! Residual vector.
            real(dp), intent(out) :: J(:,:) !! Residual Jacobian.
         end subroutine fun
      end interface

      ! ------ Arguments ------
      real(dp), intent(inout) :: x(:) !! Initial guess; last iterate on return.
      real(dp), intent(in)    :: tol !! Absolute residual convergence tolerance.
      integer, intent(in)     :: max_iter !! Maximum Newton iterations.
      integer, intent(out)    :: its !! Iteration count at termination.
      integer, intent(out)    :: info !! 0 = converged (max-norm of residuals below `tol`), 1 = iteration limit reached

      ! ------ Local variables ------
      integer :: n !! Number of unknowns
      integer :: k !! Loop/index or term index.
      real(dp) :: f_norm
      real(dp) :: phi_x
      real(dp) :: phi_trial
      real(dp) :: grad_phi_p
      real(dp) :: alpha !! Backtracking step length
      real(dp) :: f(size(X)) !! Residual or scalar equation value.
      real(dp) :: J(size(X),size(X)) !! Jacobian matrix
      real(dp) :: dx(size(X))
      real(dp) :: xtrial(size(X))
      real(dp) :: ftrial(size(X))
      real(dp) :: Jp(size(X))

      n = size(x)

      ! Assume failure until success
      info = 1

      ! ---------- Main Newton iteration ----------
      do its = 1, max_iter

         ! Evaluate f and J at current point
         call fun(x, f, J)
         f_norm = maxval(abs(f))

         if (f_norm < tol) then
            info = 0
            return
         end if

         ! Solve: J dx = -f   (replace with your solver)
         dx = solve_system(J, -f)

         ! Merit function φ = 1/2 ||f||^2
         phi_x = 0.5_dp * dot_product(f, f)

         ! Directional derivative φ' = fᵀ (J dx)
         Jp = matmul(J, dx)
         grad_phi_p = dot_product(f, Jp)

         ! If not descent, flip direction
         if (grad_phi_p > 0._dp) then
            dx = -dx
            Jp = -Jp
            grad_phi_p = -grad_phi_p
         end if

         ! ------------ Armijo backtracking ------------
         alpha = 1.0_dp

         do
            xtrial = x + alpha * dx
            call fun(xtrial, ftrial, J) ! J not needed but required by interface

            phi_trial = 0.5_dp * dot_product(ftrial, ftrial)

            if (phi_trial <= phi_x + 1e-4_dp * alpha * grad_phi_p) exit

            alpha = alpha * 0.5_dp
            if (alpha < 1e-8_dp) exit
         end do

         ! Update
         x = x + alpha * dx
      end do
   end subroutine newton_solve
end module nonlinear_solvers

module yaeos__m_s_sp
   !! # Saturation-point nonlinear systems
   !! Residual equations, Jacobians, and point solvers for bubble, dew, and liquid-liquid states.
   !!
   !! # Description
   !! The routines in this module encode saturation unknowns in logarithmic
   !! coordinates and enforce isofugacity, material-balance, and specification
   !! equations. They are implementation-level solvers used by the public
   !! `saturation_pressure` and `saturation_temperature` routines.
   !!
   !! # Examples
   !! Application code should call the public saturation-point functions rather
   !! than these solver helpers.
   use yaeos__constants, only: pr
   use yaeos__models_ar, only: ArModel, size
   implicit none

contains

   subroutine saturation_F(model, z, X, ns, S, F, dF, dPdVz, dPdVy)
      !! # saturation_F
      !! Evaluate the saturation residual system in log-volume coordinates.
      !!
      !! # Description
      !! With \(X = [\ln K_i, \ln V_z, \ln V_y, \ln T, \ln P]\) and
      !! \(y = K z\), the residuals are
      !! \(\ln f_i(y, V_y, T) - \ln f_i(z, V_z, T)\) (isofugacity),
      !! \(\sum_i (y_i - z_i)\), \(P_y - P_z\), \(\ln P - \ln P_y\), and
      !! \(X_{ns} - S\). The Jacobian is returned when `dF` is present.
      class(ArModel), intent(in) :: model !! EOS used to evaluate phase fugacities.
      real(pr), intent(in) :: X(:) !! Unknowns \([\ln K_i, \ln V_z, \ln V_y, \ln T, \ln P]\)
      integer, intent(in) :: ns !! Index of the specified unknown.
      real(pr), intent(in) :: S !! Value imposed on unknown `ns`.
      real(pr), intent(out) :: F(:) !! Saturation residual vector.
      real(pr), optional, intent(out) :: dF(:, :) !! Residual Jacobian.
      real(pr), intent(out) :: dPdVz !! Pressure derivative with respect to the feed-phase volume.
      real(pr), intent(out) :: dPdVy !! Pressure derivative with respect to the incipient-phase volume.

      ! Variables
      real(pr) :: T
      real(pr) :: Vz
      real(pr) :: Vy
      real(pr) :: z(size(model)) !! Feed composition

      ! Main phase variables
      real(pr) :: lnfug_z(size(model)) !! Logarithmic thermodynamic quantity or intermediate.
      real(pr) :: dlnfug_dn_z(size(model), size(model))
      real(pr) :: dlnfug_dT_z(size(model))
      real(pr) :: dlnfug_dV_z(size(model))
      real(pr) :: dlnfug_dP_z(size(model))
      real(pr) :: Pz
      real(pr) :: dPdTz
      real(pr) :: dPdn_z(size(z))

      ! incipient phase variables
      real(pr) :: y(size(z)) !! Incipient/second phase composition.
      real(pr) :: lnfug_y(size(model)) !! Logarithmic thermodynamic quantity or intermediate.
      real(pr) :: dlnfug_dn_y(size(model), size(model))
      real(pr) :: dlnfug_dT_y(size(model))
      real(pr) :: dlnfug_dV_y(size(model))
      real(pr) :: dlnfug_dP_y(size(model))
      real(pr) :: Py
      real(pr) :: dPdTy
      real(pr) :: dPdn_y(size(z))

      real(pr) :: lnPspec !! Logarithmic thermodynamic quantity or intermediate.

      integer :: j !! Loop/index variable.
      integer :: nc !! Number of components.

      nc = size(z)

      y  = z * exp(X(:nc))
      Vz = exp(X(nc+1))
      Vy = exp(X(nc+2))
      T  = exp(X(nc+3))
      lnPspec = X(nc+4)

      if (present(df)) then
         call model%lnfug_vt(&
            n=z, V=Vz, T=T, P=Pz, dPdT=dPdTz, dPdV=dPdVz, dPdn=dPdn_z, &
            lnf=lnfug_z, &
            dlnfdV=dlnfug_dV_z, dlnfdT=dlnfug_dT_z, dlnfdn=dlnfug_dn_z &
            )

         call model%lnfug_vt(&
            n=y, V=Vy, T=T, P=Py, dPdT=dPdTy, dPdV=dPdVy, dPdn=dPdn_y, &
            lnf=lnfug_y, &
            dlnfdV=dlnfug_dV_y, dlnfdT=dlnfug_dT_y, dlnfdn=dlnfug_dn_y &
            )
      else
         call model%lnfug_vt(n=z, V=Vz, T=T, P=Pz, lnf=lnfug_z)
         call model%lnfug_vt(n=y, V=Vy, T=T, P=Py, lnf=lnfug_y)
      end if

      F = 0

      F(:nc) = lnfug_y - lnfug_z
      F(nc + 1) = sum(y - z)
      F(nc + 2) = Py - Pz
      F(nc + 3) = lnPspec - log(Py)
      F(nc + 4) = X(ns) - S

      if (present(dF)) then
         dF = 0

         ! isofugacity
         do j=1,nc
            df(:nc, j) = y(j) * dlnfug_dn_y(:, j)
         end do

         dF(:nc, nc+1) = -dlnfug_dV_z * Vz
         dF(:nc, nc+2) =  dlnfug_dV_y * Vy
         dF(:nc, nc+3) =  T * (dlnfug_dT_y - dlnfug_dT_z)

         ! mass balance
         df(nc+1, :nc) = y

         ! pressure equality
         df(nc+2, :nc)  = y * dPdn_y
         df(nc+2, nc+1) = -dPdVz * Vz
         df(nc+2, nc+2) = dPdVy * Vy
         df(nc+2, nc+3) = T*(dPdTy - dPdTz)

         df(nc+3, :nc) = -y * dPdn_y/Py
         df(nc+3, nc+1) = 0
         df(nc+3, nc+2) = -dPdVy * Vy / Py
         df(nc+3, nc+3) = -dPdTy * T / Py
         df(nc+3, nc+4) = 1

         df(nc+4, ns) = 1
      end if
   end subroutine saturation_F

   subroutine saturation_TP(model, kind, z, X, ns, S, F, dF, dFdS)
      !! # saturation_TP
      !! Evaluate a saturation residual system at fixed pressure or temperature.
      !!
      !! # Description
      !! Reconstructs phase compositions from logarithmic K-factors, evaluates
      !! fugacity coefficients with root selections determined by `kind`, and
      !! returns isofugacity, material-balance, and specification residuals:
      !! \(\ln K_i + \ln\phi_i(y) - \ln\phi_i(z)\), \(\sum_i (y_i - z_i)\),
      !! and \(X_{ns} - S\).
      class(ArModel), intent(in) :: model !! EOS used for phase properties.
      character(len=*), intent(in) :: kind !! `"bubble"`, `"dew"`, or `"liquid-liquid"`; any other value uses stable roots.
      real(pr), intent(in) :: z(size(model)) !! Feed composition.

      real(pr), intent(in) :: X(:) !! Unknowns \([\ln K_i, \ln T, \ln P]\)
      integer, intent(in) :: ns !! Index of the specified unknown.
      real(pr), intent(in) :: S !! Value imposed on unknown `ns`.

      real(pr), intent(out) :: F(:) !! Residual vector.
      real(pr), intent(out) :: dF(:, :) !! Jacobian with respect to the unknown vector.
      real(pr), intent(out) :: dFdS(:) !! Residual derivative with respect to the specification.

      character(len=14) :: kind_z
      character(len=14) :: kind_y

      real(pr) :: y(size(X)-2) !! Incipient/second phase composition.
      real(pr) :: Vz
      real(pr) :: Vy
      real(pr) :: lnPhi_z(size(X)-2) !! Logarithmic thermodynamic quantity or intermediate.
      real(pr) :: lnPhi_y(size(X)-2) !! Logarithmic thermodynamic quantity or intermediate.
      real(pr) :: dlnphi_dt_z(size(X)-2)
      real(pr) :: dlnphi_dt_y(size(X)-2)
      real(pr) :: dlnphi_dp_z(size(X)-2)
      real(pr) :: dlnphi_dp_y(size(X)-2)
      real(pr) :: dlnphi_dn_z(size(X)-2, size(X)-2)
      real(pr) :: dlnphi_dn_y(size(X)-2, size(X)-2)

      real(pr) :: T
      real(pr) :: P
      real(pr) :: K(size(X)-2) !! Equilibrium K-factors

      integer :: i !! Loop index `i`.
      integer :: j !! Loop index `j`.
      integer :: nc !! Number of components.

      nc = size(X)-2

      F = 0
      dF = 0

      K = exp(X(:nc))
      T = exp(X(nc+1))
      P = exp(X(nc+2))

      y = K*z

      select case(kind)
       case ("bubble")
         kind_z = "liquid"
         kind_y = "vapor"
       case ("dew")
         kind_z = "vapor"
         kind_y = "liquid"
       case ("liquid-liquid")
         kind_z = "liquid"
         kind_y = "liquid"
       case default
         kind_z = "stable"
         kind_y = "stable"
      end select

      call model%lnphi_pt(&
         z, P, T, V=Vz, root_type=kind_z, &
         lnPhi=lnphi_z, dlnPhidt=dlnphi_dt_z, &
         dlnPhidp=dlnphi_dp_z, dlnphidn=dlnphi_dn_z &
         )
      call model%lnphi_pt(&
         y, P, T, V=Vy, root_type=kind_y, &
         lnPhi=lnphi_y, dlnPhidt=dlnphi_dt_y, &
         dlnPhidp=dlnphi_dp_y, dlnphidn=dlnphi_dn_y &
         )

      F(:nc) = X(:nc) + lnPhi_y - lnPhi_z
      F(nc + 1) = sum(y - z)
      F(nc + 2) = X(ns) - S

      ! Jacobian Matrix
      do j=1,nc
         df(:nc, j) = dlnphi_dn_y(:, j) * y(j)
         df(j, j) = dF(j, j) + 1
      end do

      df(:nc, nc + 1) = T * (dlnphi_dt_y - dlnphi_dt_z)
      df(:nc, nc + 2) = P * (dlnphi_dp_y - dlnphi_dp_z)

      df(nc + 1, :nc) = y

      df(nc + 2, :) = 0
      df(nc + 2, ns) = 1

      dFdS = 0
      dFdS(nc+2) = -1
   end subroutine saturation_TP

   subroutine solve_TP(model, kind, z, X, ns, S, tol, max_iterations, its)
      !! # solve_TP
      !! Correct a fixed-pressure or fixed-temperature saturation point.
      !!
      !! # Description
      !! Applies `newton` from `yaeos__math_nonlinearsolvers` to the residual
      !! system assembled by `saturation_TP` and updates the packed unknown
      !! vector in place.
      use nonlinear_solvers, only: newton_solve
      use yaeos__math, only: solve_system
      use yaeos__math_nonlinearsolvers, only: newton, homotopy
      class(ArModel), intent(in) :: model !! EOS used for phase fugacities.
      character(len=*), intent(in) :: kind !! `"bubble"`, `"dew"`, or `"liquid-liquid"`; any other value uses stable roots.
      real(pr), intent(in) :: z(:) !! Feed composition.
      real(pr), intent(inout) :: X(:) !! Initial packed unknowns; corrected state on return.
      integer, intent(in) :: ns !! Index of the specified variable.
      real(pr), intent(in) :: S !! Specified value of `X(ns)` (log coordinates)
      real(pr), intent(in) :: tol !! Residual tolerance.
      integer, intent(in) :: max_iterations !! Maximum corrector iterations.
      integer, intent(out) :: its !! Iterations used.

      integer :: nc !! Number of components.

      real(pr) :: F(size(X)) !! Residual or scalar equation value.
      real(pr) :: dF(size(X), size(X))
      real(pr) :: dFdS(size(X))
      real(pr) :: dx(size(X))

      real(pr) :: t
      real(pr) :: G(size(X))
      real(pr) :: dG(size(X), size(X))
      real(pr) :: H(size(X))
      real(pr) :: dH(size(X), size(X))
      real(pr) :: X0(size(X))

      real(pr) :: alpha !! Step length (unused)
      real(pr) :: phi_x
      real(pr) :: phi_ax
      real(pr) :: grad_phi_p
      real(pr) :: Jp(size(X))
      real(pr) :: Xnew(size(X))
      integer :: info !! Solver status code.

      nc = size(X) - 2

      ! call homotopy(sub=wrap, x=X, tol=tol, max_its=max_iterations, its=its)
      call newton(sub=wrap, x=X, tol=tol, max_its=max_iterations, its=its)

      ! call newton_solve(fun=wrap, x=X, tol=1e-7_pr, max_iter=20, its=its, info=info)
      call wrap(X, F, dF)
      ! print *, F
      ! print *, its, maxval(abs(F))
   contains
      subroutine wrap(X, F, J)
         real(pr), intent(in) :: X(:) !! Unknowns \([\ln K_i, \ln T, \ln P]\)
         real(pr), intent(out) :: F(:) !! Residual vector
         real(pr), intent(out) :: J(:, :) !! Jacobian matrix
         call saturation_TP(model, kind, z, X, ns, S=S, F=F, dF=J, dFdS=dFdS)
      end subroutine
   end subroutine solve_TP

   subroutine solve_VxVyT(model, z, X, ns, S, tol, max_iterations, its)
      !! # solve_VxVyT
      !! Solve saturation equations using log phase volumes and temperature.
      !!
      !! # Description
      !! Undamped Newton correction of the `saturation_F` system, with unknowns
      !! \([\ln K_i, \ln V_z, \ln V_y, \ln T, \ln P]\). Pressure is an unknown
      !! tied to the incipient-phase pressure, and phase pressure equality is
      !! one of the residuals.
      use yaeos__math, only: solve_system
      class(ArModel), intent(in) :: model !! EOS used for phase properties.
      real(pr), intent(in) :: z(:) !! Feed composition.
      real(pr), intent(inout) :: X(:) !! Initial unknown vector; corrected state on return.
      integer, intent(in) :: ns !! Index of the specified variable.
      real(pr), intent(in) :: S !! Specification value in solver coordinates.
      real(pr), intent(in) :: tol !! Residual tolerance.
      integer, intent(in) :: max_iterations !! Maximum Newton iterations.
      integer, intent(out) :: its !! Iterations used.

      real(pr) :: dPdVz !! Pressure derivative with respect to the feed-phase volume.
      real(pr) :: dPdVy !! Pressure derivative with respect to the incipient-phase volume.

      integer :: nc !! Number of components.
      real(pr) :: F(size(X)) !! Residual or scalar equation value.
      real(pr) :: dF(size(X), size(X))
      real(pr) :: dFdS(size(X))
      real(pr) :: Xold(size(X))
      real(pr) :: dx(size(X))
      real(pr) :: dx_old(size(x))

      nc = size(X) - 4

      its = 0
      do while (its < max_iterations)
         call saturation_F(model, z, X, ns, S, F, dF, dPdVz, dPdVy)
         if (all(abs(F) < tol)) exit

         dX = solve_system(dF, -F)

         X = X + dX

         its = its + 1
      end do
   end subroutine solve_VxVyT
end module yaeos__m_s_sp
