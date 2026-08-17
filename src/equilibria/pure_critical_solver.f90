module critical_point_solver
   !! # `critical_point_solver`
   !! Calculation of pure-component critical points 
   !! using the yaeos equation of state and numerical derivatives.
   !!
   !! # Description
   !! The critical point of a pure component is defined by the conditions:
   !!
   !! \[
   !! \frac{\partial P}{\partial V} = 0,
   !! \frac{\partial^2 P}{\partial V^2} = 0
   !! \]
   !!
   !! This module solves the resulting nonlinear system in the variables
   !! `(V, T)` using a Newton-Raphson iteration. Each residual is evaluated
   !! from the equation-of-state pressure and its numerical derivatives.
   !!
   !! # Notes
   !! The implementation is intentionally simple and robust for pure-component
   !! searches, using finite differences and a safeguarded Newton step.
   use yaeos__constants, only: pr
   use yaeos__models_ar, only: ArModel

   implicit none
   private
   public :: find_critical_points_all_components

   integer, parameter :: MAX_ITER = 100
   real(pr), parameter :: TOL_RES = 1.0e-8_pr     ! residual tolerance
   real(pr), parameter :: TOL_STEP = 1.0e-6_pr    ! step tolerance

   ! Steps for numerical derivatives
   real(pr), parameter :: DV_STEP = 1.0e-4_pr     ! L/mol
   real(pr), parameter :: DT_STEP = 1.0e-3_pr     ! K

contains

   ! ========================================================================
   ! MAIN ROUTINE: Find critical point for each pure component
   ! ========================================================================

   subroutine find_critical_points_all_components(model, nc, Vc, Tc, Pc, converged)
      !! # `find_critical_points_all_components`
      !! Calculate the critical point of each pure component in a mixture.
      !!
      !! # Description
      !! This routine estimates an initial point for each component and solves the
      !! criticality equations for the corresponding pure fluid. The routine
      !! stores the resulting critical volume, temperature and pressure in the
      !! output arrays.
      !!
      !! # Arguments
      !! - `model`: Equation of state model used to evaluate pressure and
      !!   derivatives.
      !! - `nc`: Number of components.
      !! - `Vc(nc)`: Critical molar volume of each component [L/mol].
      !! - `Tc(nc)`: Critical temperature of each component [K].
      !! - `Pc(nc)`: Critical pressure of each component [bar].
      !! - `converged(nc)`: Logical flag indicating whether the iteration
      !!   converged for each component.

      class(ArModel), intent(in) :: model !! Equation of state model used to evaluate pressure and derivatives.
      integer, intent(in) :: nc !! Number of components.
      real(pr), intent(out) :: Vc(nc) !! Critical molar volume [L/mol]
      real(pr), intent(out) :: Tc(nc) !! Critical temperature [K]
      real(pr), intent(out) :: Pc(nc) !! Critical pressure [bar] for each component.
      logical, intent(out) :: converged(nc) !! Convergence flag for each component.

      integer :: i               ! component loop counter
      real(pr) :: V0, T0         ! initial critical-volume and temperature guesses
      logical :: conv            ! convergence flag for this component

      do i = 1, nc
         ! Initial estimation: based on simple correlations
         call estimate_critical_point(i, nc, V0, T0)

         ! Solve Newton-Raphson 2D
         call solve_critical_point(model, nc, i, V0, T0, Vc(i), Tc(i), Pc(i), conv)

         converged(i) = conv

      end do

   end subroutine find_critical_points_all_components


   ! ========================================================================
   ! INITIAL ESTIMATE of critical point
   ! ========================================================================

   subroutine estimate_critical_point(i, nc, V_est, T_est)
      !! # `estimate_critical_point`
      !! Produce a simple initial guess for the critical volume and temperature.
      !!
      !! # Description
      !! This helper builds a rough estimate for component `i` based on a
      !! normalized index ranging from the lightest to heaviest component. The
      !! values are intentionally heuristic and only intended to initialize the
      !! Newton solver.
      !!
      !! # Arguments
      !! - `i`: Component index.
      !! - `nc`: Total number of components.
      !! - `V_est`: Estimated critical molar volume [L/mol].
      !! - `T_est`: Estimated critical temperature [K].

      integer, intent(in) :: i !! Component index.
      integer, intent(in) :: nc !! Total number of components.
      real(pr), intent(out) :: V_est !! Estimated critical molar volume [l/mol].
      real(pr), intent(out) :: T_est !! Estimated critical temperature [K].

      real(pr) :: mw              ! molecular weight
      real(pr) :: f               ! normalized component index

      ! Normalized index: 0 (light) -> 1 (heavy)
      f = real(i - 1, pr) / max(nc - 1, 1)

      ! Molecular weight: C1 (~16) to pseudo-component (~416)
      mw = 16.0_pr + f * 400.0_pr

      ! Empirical heuristics for critical point
      ! Based on typical behavior of hydrocarbons
      ! V_c ≈ 0.27 * R * T_c / P_c (law of corresponding states)

      T_est = 300! 150.0_pr + f * 350.0_pr        ! Tc: 150 K (CH4) to 500 K (heavy)
      V_est = 1!1e-3_pr * (50.0_pr + f * 150.0_pr)         ! Vc: 50 to 200 L/mol

   end subroutine estimate_critical_point


   ! ========================================================================
   ! SOLVER: Newton-Raphson 2D to solve [dPdV, d²PdV²] = 0
   ! ========================================================================

   subroutine solve_critical_point(model, nc, i_comp, V0, T0, Vc, Tc, Pc, converged)
      !! # `solve_critical_point`
      !! Solve the pure-component criticality conditions using a Newton method.
      !!
      !! # Description
      !! The system is defined by the residual vector:
      !!
      !! \[
      !! f_1(V,T) = \frac{\partial P}{\partial V}, \quad
      !! f_2(V,T) = \frac{\partial^2 P}{\partial V^2}
      !! \]
      !!
      !! with unknowns `(V, T)`. The Jacobian is approximated numerically and a
      !! safeguarded step is applied to keep the iterates inside a physically
      !! meaningful region.
      !!
      !! # Arguments
      !! - `model`: Equation of state model.
      !! - `nc`: Number of components.
      !! - `i_comp`: Component being solved.
      !! - `V0`: Initial guess for the molar volume [L/mol].
      !! - `T0`: Initial guess for the temperature [K].
      !! - `Vc`: Final critical volume [L/mol].
      !! - `Tc`: Final critical temperature [K].
      !! - `Pc`: Final critical pressure [bar].
      !! - `converged`: `.true.` if the critical point was found.

      class(ArModel), intent(in) :: model !! Equation of state model.
      integer, intent(in) :: nc !! Number of components.
      integer, intent(in) :: i_comp !! Component being solved.
      real(pr), intent(in) :: V0 !! Initial guess for the molar volume [L/mol].
      real(pr), intent(in) :: T0 !! Initial guess for the temperature [K].
      real(pr), intent(out) :: Vc !! Final critical volume [L/mol].
      real(pr), intent(out) :: Tc !! Final critical temperature [K].
      real(pr), intent(out) :: Pc !! Final critical pressure [bar].
      logical, intent(out) :: converged !! `.true.` if the critical point was found.

      real(pr) :: V, T, P, dV, dT ! current volume, temperature, pressure, and step sizes
      real(pr) :: f1, f2, norm_f ! residuals and residual norm
      real(pr) :: J(2,2)         ! Jacobian matrix
      real(pr) :: J_inv(2,2)     ! inverse Jacobian matrix
      real(pr) :: det            ! determinant
      real(pr) :: dX(2)          ! Newton step
      integer :: iter            ! Newton iteration counter

      V = V0
      T = T0
      converged = .false.

      ! Newton-Raphson loop
      do iter = 1, MAX_ITER

         ! Evaluate residuals: f1 = dP/dV,  f2 = d²P/dV²
         call residual_critical_eqs(model, nc, i_comp, V, T, f1, f2, P)

         norm_f = sqrt(f1**2 + f2**2)

         ! Convergence criterion
         if (norm_f < TOL_RES) then
            converged = .true.
            Vc = V
            Tc = T
            Pc = P
            return
         end if

         ! Calculate Jacobian numerically
         call jacobian_critical_eqs(model, nc, i_comp, V, T, J)

         ! Invert 2x2 matrix
         det = J(1,1) * J(2,2) - J(1,2) * J(2,1)
         if (abs(det) < 1.0e-12_pr) then
            ! Singular Jacobian, abort
            Vc = V
            Tc = T
            return
         end if

         J_inv(1,1) = J(2,2) / det
         J_inv(1,2) = -J(1,2) / det
         J_inv(2,1) = -J(2,1) / det
         J_inv(2,2) = J(1,1) / det

         ! Newton step: dX = -J_inv * f
         dX(1) = -(J_inv(1,1) * f1 + J_inv(1,2) * f2)
         dX(2) = -(J_inv(2,1) * f1 + J_inv(2,2) * f2)

         ! Update variables
         V = V + dX(1)
         T = T + dX(2)

         ! Safeguards: keep V, T in physical ranges
         if (V <= 0.0_pr) V = V0 / 2.0_pr   ! bounce if V < 0
         if (T <= 0.0_pr) T = T0 / 2.0_pr   ! bounce if T < 0

         ! Step convergence criterion
         if (sqrt(dX(1)**2 + dX(2)**2) < TOL_STEP .and. norm_f < 10.0_pr * TOL_RES) then
            converged = .true.
            Vc = V
            Tc = T
            return
         end if

      end do

      ! If no convergence after MAX_ITER
      Vc = V
      Tc = T
      converged = .false.

   end subroutine solve_critical_point


   ! ========================================================================
   ! Evaluate residuals: [f1, f2] = [dP/dV, d²P/dV²]
   ! ========================================================================

   subroutine residual_critical_eqs(model, nc, i_comp, V, T, f1, f2, P)
      !! # `residual_critical_eqs`
      !! Evaluate the residuals of the criticality equations.
      !!
      !! # Description
      !! This routine evaluates the pressure and the first derivative with the
      !! EOS, then approximates the second derivative with respect to volume by a
      !! centered finite difference. The resulting residuals correspond to the
      !! conditions at the critical point.
      !!
      !! # Arguments
      !! - `model`: Equation of state model.
      !! - `nc`: Number of components.
      !! - `i_comp`: Component index.
      !! - `V`: Molar volume [L/mol].
      !! - `T`: Temperature [K].
      !! - `f1`: Residual associated with `dP/dV`.
      !! - `f2`: Residual associated with `d²P/dV²`.
      !! - `P`: Pressure evaluated at `(V,T)` [bar].

      class(ArModel), intent(in) :: model !! Equation of state model.
      integer, intent(in) :: i_comp !! Component index.
      integer, intent(in) :: nc !! Number of components.
      real(pr), intent(in) :: V !! Molar volume [L/mol].
      real(pr), intent(in) :: T !! Temperature [K].
      real(pr), intent(out) :: f1 !! Residual for `dP/dV`.
      real(pr), intent(out) :: f2 !! Residual for `d²P/dV²`.
      real(pr), intent(out) :: P  !! Pressure [bar].

      real(pr) :: z(nc)            ! composition vector
      real(pr) :: dPdV             ! pressure derivative
      real(pr) :: dPdV_plus        ! forward finite-difference derivative
      real(pr) :: dPdV_minus       ! backward finite-difference derivative

      ! Pure composition for component i_comp
      z = 0.0_pr
      z(i_comp) = 1.0_pr

      ! First derivative from yaeos
      call model%pressure(z, V, T, P, dPdV=dPdV)

      f1 = dPdV

      ! Second derivative via centered numerical difference
      ! d²P/dV² ≈ [dP/dV(V+h) - dP/dV(V-h)] / (2*h)
      call model%pressure(z, V + DV_STEP, T, P, dPdV=dPdV_plus)
      call model%pressure(z, V - DV_STEP, T, P, dPdV=dPdV_minus)

      f2 = (dPdV_plus - dPdV_minus) / (2.0_pr * DV_STEP)

   end subroutine residual_critical_eqs


   ! ========================================================================
   ! Jacobian: matrix of derivatives ∂f_i / ∂x_j
   ! ========================================================================

   subroutine jacobian_critical_eqs(model, nc, i_comp, V, T, J)
      !! # `jacobian_critical_eqs`
      !! Approximate the Jacobian of the criticality residuals.
      !!
      !! # Description
      !! The Jacobian is computed by centered finite differences:
      !!
      !! \[
      !! J = \begin{bmatrix}
      !! \partial f_1/\partial V & \partial f_1/\partial T \\
      !! \partial f_2/\partial V & \partial f_2/\partial T
      !! \end{bmatrix}
      !! \]
      !!
      !! where `f1 = dP/dV` and `f2 = d²P/dV²`. This numerical approximation is
      !! used in the Newton update for `(V, T)`.
      !!
      !! # Arguments
      !! - `model`: Equation of state model.
      !! - `nc`: Number of components.
      !! - `i_comp`: Component index.
      !! - `V`: Molar volume [L/mol].
      !! - `T`: Temperature [K].
      !! - `J(2,2)`: Numerical Jacobian matrix.

      class(ArModel), intent(in) :: model !! Equation of state model.
      integer, intent(in) :: i_comp !! Component index.
      integer, intent(in) :: nc !! Number of components.
      real(pr), intent(in) :: V !! Molar volume [L/mol].
      real(pr), intent(in) :: T !! Temperature [K].
      real(pr), intent(out) :: J(2,2) !! Numerical Jacobian matrix of the criticality residuals.

      real(pr) :: f1_V_plus  ! f1 finite-difference evaluation at V+h
      real(pr) :: f1_V_minus ! f1 finite-difference evaluation at V-h
      real(pr) :: f1_T_plus  ! f1 finite-difference evaluation at T+h
      real(pr) :: f1_T_minus ! f1 finite-difference evaluation at T-h
      real(pr) :: f2_V_plus  ! f2 finite-difference evaluation at V+h
      real(pr) :: f2_V_minus ! f2 finite-difference evaluation at V-h
      real(pr) :: f2_T_plus  ! f2 finite-difference evaluation at T+h
      real(pr) :: f2_T_minus ! f2 finite-difference evaluation at T-h

      real(pr) :: P_V_plus  ! pressure at V+h
      real(pr) :: P_V_minus ! pressure at V-h
      real(pr) :: P_T_plus  ! pressure at T+h
      real(pr) :: P_T_minus ! pressure at T-h

      ! J(1,1) = ∂f1/∂V ≈ [f1(V+h) - f1(V-h)] / (2*h)
      call residual_critical_eqs(model, nc, i_comp, V + DV_STEP, T, f1_V_plus, f2_V_plus, P_V_plus)
      call residual_critical_eqs(model, nc, i_comp, V - DV_STEP, T, f1_V_minus, f2_V_minus, P_V_minus)
      J(1,1) = (f1_V_plus - f1_V_minus) / (2.0_pr * DV_STEP)

      ! J(1,2) = ∂f1/∂T ≈ [f1(T+h) - f1(T-h)] / (2*h)
      call residual_critical_eqs(model, nc, i_comp, V, T + DT_STEP, f1_T_plus, f2_T_plus, P_T_plus)
      call residual_critical_eqs(model, nc, i_comp, V, T - DT_STEP, f1_T_minus, f2_T_minus, P_T_minus)
      J(1,2) = (f1_T_plus - f1_T_minus) / (2.0_pr * DT_STEP)

      ! J(2,1) = ∂f2/∂V ≈ [f2(V+h) - f2(V-h)] / (2*h)
      J(2,1) = (f2_V_plus - f2_V_minus) / (2.0_pr * DV_STEP)

      ! J(2,2) = ∂f2/∂T ≈ [f2(T+h) - f2(T-h)] / (2*h)
      J(2,2) = (f2_T_plus - f2_T_minus) / (2.0_pr * DT_STEP)

   end subroutine jacobian_critical_eqs

end module critical_point_solver
