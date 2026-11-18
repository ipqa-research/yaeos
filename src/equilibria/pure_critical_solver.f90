module yaeos__critical_pure_point_solver
   !! # Pure-component critical-point solver
   !! Find pure-fluid critical constants for each component in an EOS model.
   !!
   !! # Description
   !! The critical point of a pure component is defined by the conditions:
   !!
   !! \[
   !! \frac{\partial P}{\partial V} = 0, \quad \frac{\partial^2 P}{\partial V^2} = 0
   !! \]
   !!
   !! This module solves the resulting nonlinear system in the variables
   !! `(V, T)` using a Newton-Raphson iteration. Each residual is evaluated
   !! from the equation-of-state pressure and its numerical derivatives.
   !!
   !! # Examples
   !! ```fortran
   !! call find_critical_points_all_components(model, nc, Vc, Tc, Pc, converged)
   !! ```

   use yaeos__constants, only: pr
   use yaeos__models_ar, only: ArModel

   implicit none
   private
   public :: find_critical_points_all_components

   integer, parameter :: MAX_ITER = 151 !! Maximum Newton iterations per component.
   real(pr), parameter :: TOL_RES = 2.0e-8_pr !! Convergence threshold for the critical residual norm.
   real(pr), parameter :: TOL_STEP = 2.0e-6_pr !! Convergence threshold for the Newton step norm.

   ! Steps for numerical derivatives
   real(pr), parameter :: DV_STEP = 2.0e-4_pr !! Finite-difference increment in volume [L/mol].
   real(pr), parameter :: DT_STEP = 2.0e-2_pr !! Finite-difference increment in temperature [K].

contains

   ! ========================================================================
   ! MAIN ROUTINE: Find critical point for each pure component
   ! ========================================================================
   subroutine find_critical_points_all_components(model, nc, Vc, Tc, Pc, converged)
      !! # find_critical_points_all_components
      !! Estimate pure-component critical constants for every model component.
      !!
      !! # Description
      !! For each component, first obtains a starting estimate from pure-fluid
      !! saturation calculations, then solves the criticality equations
      !! \(P_V=0\) and \(P_{VV}=0\) with a two-variable Newton iteration.
      !! A failed solve is reported in `converged`; the corresponding output
      !! values contain the final iterate and must not be treated as validated
      !! critical constants.
      !!
      !! # Examples
      !! ```fortran
      !! call find_critical_points_all_components(model, nc, Vc, Tc, Pc, ok)
      !! ```
      class(ArModel), intent(in) :: model !! Equation of state used for pure-fluid saturation properties.
      integer, intent(in) :: nc !! Number of components; must match the model.
      real(pr), intent(out) :: Vc(nc) !! Critical molar volumes [L/mol].
      real(pr), intent(out) :: Tc(nc) !! Critical temperatures [K].
      real(pr), intent(out) :: Pc(nc) !! Critical pressures [bar].
      logical, intent(out) :: converged(nc) !! True for components whose solve converged.

      integer :: i !! Loop/index variable.
      real(pr) :: V_init
      real(pr) :: T_init
      logical :: conv

      ! Fixed: Loop should usually start at 1 to cover all components
      do i = 1, nc
         call estimate_critical_point(model, i, nc, V_init, T_init)
         call solve_critical_point(model, nc, i, V_init, T_init, Vc(i), Tc(i), Pc(i), conv)
         converged(i) = conv
      end do
   end subroutine find_critical_points_all_components


   ! ========================================================================
   ! INITIAL ESTIMATE of critical point
   ! ========================================================================
   subroutine estimate_critical_point(model, i, nc, V_est, T_est)
      !! # estimate_critical_point
      !! Construct an initial pure-component critical-point estimate.
      !!
      !! # Description
      !! Uses successive pure-fluid saturation calculations and the ratio of
      !! vapor to liquid volume to move the trial temperature toward the
      !! critical region. The estimates are initial guesses, not final results.
      use yaeos__equilibria_boundaries_pure_saturation, only:pure_saturation_line, PurePsat
      class(ArModel), intent(in) :: model !! Equation of state used to evaluate pressure and derivatives.
      integer, intent(in) :: i !! One-based component index.
      integer, intent(in) :: nc !! Number of components in the model.
      real(pr), intent(out) :: V_est !! Estimated critical volume [L/mol].
      real(pr), intent(out) :: T_est !! Estimated critical temperature [K].

      real(pr) :: z(nc) !! Feed or overall composition.
      real(pr) :: T
      real(pr) :: P_sat
      real(pr) :: VL
      real(pr) :: VV
      real(pr) :: ratio
      real(pr) :: dT
      logical :: psat_converged
      integer :: iter !! Iteration count.

      type(PurePsat) :: psat !! Declared saturation-line result; currently unused by this estimate.

      z = 0.0_pr
      z(i) = 1.0_pr

      ! Initial temperature heuristic scan
      T = 150.0_pr
      V_est = 0.15_pr
      T_est = T


      iter = 0
      do while(.true. .and. iter < 100)
         P_sat = model%Psat_pure(i, T, Vl=VL, Vv=Vv, converged=psat_converged)
         iter = iter + 1
         if (abs(Vv - Vl) > 1e-3) then
            ratio = Vv / VL
            if (ratio < 10_pr) then
               T_est = T
               V_est = Vl
               return
            end if
            ! Advance temperature upward towards Tc
            dT = 0.08_pr * log(ratio)
            dT = sign(min(abs(dT), 0.1_pr), dT)
            T = T  * (1.0_pr + dT)
         else
            ! Step temperature down if above pseudo-critical or single phase
            T = T * 0.75
         end if
      end do
   end subroutine estimate_critical_point


   ! ========================================================================
   ! SOLVER: Newton-Raphson 2D to solve [dPdV, d²PdV²] = 0
   ! ========================================================================
   subroutine solve_critical_point(model, nc, i_comp, V_init, T_init, Vc, Tc, Pc, converged)
      !! # solve_critical_point
      !! Solve the pure-component criticality equations from an initial guess.
      !!
      !! # Description
      !! The residual vector is \((P_V,P_{VV})\). Its Jacobian with respect to
      !! \((V,T)\) is approximated by centered finite differences, and a direct
      !! 2-by-2 Newton step is applied. Non-positive trial volume or temperature
      !! is reset to half the corresponding initial value.
      class(ArModel), intent(in) :: model !! Equation of state used for pressure derivatives.
      integer, intent(in) :: nc !! Number of components.
      integer, intent(in) :: i_comp !! One-based index of the pure component.
      real(pr), intent(in) :: V_init !! Initial molar volume [L/mol].
      real(pr), intent(in) :: T_init !! Initial temperature [K].
      real(pr), intent(out) :: Vc !! Critical molar volume [L/mol].
      real(pr), intent(out) :: Tc !! Critical temperature [K].
      real(pr), intent(out) :: Pc !! Critical pressure [bar].
      logical, intent(out) :: converged !! Whether either convergence test succeeded.

      real(pr) :: V
      real(pr) :: T
      real(pr) :: P
      real(pr) :: f1
      real(pr) :: f2
      real(pr) :: norm_f
      real(pr) :: J(2,2) !! Jacobian matrix.
      real(pr) :: J_inv(2,2)
      real(pr) :: det
      real(pr) :: dX(2)
      integer :: iter !! Iteration count.

      V = V_init
      T = T_init
      converged = .false.

      ! Newton-Raphson loop
      do iter = 1, MAX_ITER
         ! Evaluate residuals: f1 = dP/dV,  f2 = d²P/dV²
         call residual_critical_eqs(model, nc, i_comp, V, T, f1, f2, P)

         norm_f = sqrt(f1**2 + f2**2)

         if (norm_f < TOL_RES) then
            converged = .true.
            Vc = V; Tc = T; Pc = P
            return
         end if

         ! Calculate 2x2 Jacobian numerically
         call jacobian_critical_eqs(model, nc, i_comp, V, T, J)

         ! Invert 2x2 matrix
         det = J(1,1) * J(2,2) - J(1,2) * J(2,1)
         if (abs(det) < 2.0e-12_pr) then
            Vc = V; Tc = T; Pc = P
            return ! Singular Jacobian
         end if

         J_inv(1,1) =  J(2,2) / det
         J_inv(1,2) = -J(1,2) / det
         J_inv(2,1) = -J(2,1) / det
         J_inv(2,2) =  J(1,1) / det

         ! Newton step: dX = -J_inv * f
         dX(1) = -(J_inv(1,1) * f1 + J_inv(1,2) * f2)
         dX(2) = -(J_inv(2,1) * f1 + J_inv(2,2) * f2)

         ! Update variables
         V = V + dX(1)
         T = T + dX(2)

         ! Safeguards
         if (V <= 0.0_pr) V = V_init / 2.0_pr
         if (T <= 0.0_pr) T = T_init / 2.0_pr

         ! Step convergence criterion
         if (sqrt(dX(1)**2 + dX(2)**2) < TOL_STEP .and. norm_f < 10.0_pr * TOL_RES) then
            converged = .true.
            Vc = V; Tc = T; Pc = P
            return
         end if
      end do

      Vc = V; Tc = T; Pc = P
      converged = .false.
   end subroutine solve_critical_point


   ! ========================================================================
   ! Evaluate residuals: [f1, f2] = [dP/dV, d²P/dV²]
   ! ========================================================================
   subroutine residual_critical_eqs(model, nc, i_comp, V, T, f1, f2, P)
      !! # residual_critical_eqs
      !! Evaluate the two criticality residuals at a pure-fluid state.
      !!
      !! # Description
      !! Computes \(f_1=P_V\) from the EOS and \(f_2=P_{VV}\) by a centered
      !! finite difference of \(P_V\) in volume.
      class(ArModel), intent(in) :: model !! Equation of state used for pressure derivatives.
      integer, intent(in) :: i_comp !! One-based index of the pure component.
      integer, intent(in) :: nc !! Number of components.
      real(pr), intent(in) :: V !! Total volume [L].
      real(pr), intent(in) :: T !! Temperature [K].
      real(pr), intent(out) :: f1 !! First criticality residual, pressure differentiated with respect to volume.
      real(pr), intent(out) :: f2 !! Second criticality residual, pressure differentiated twice with respect to volume.
      real(pr), intent(out) :: P !! Pressure [bar].

      real(pr) :: z(nc)
      real(pr) :: dPdV
      real(pr) :: dPdV_plus
      real(pr) :: dPdV_minus
      real(pr) :: Pin

      ! Fixed: Pure component mole fraction is 1.0, others 0.0
      z = 0.0_pr
      z(i_comp) = 1.0_pr

      call model%pressure(z, V, T, P, dPdV=dPdV)
      f1 = dPdV

      call model%pressure(z, V + DV_STEP, T, Pin, dPdV=dPdV_plus)
      call model%pressure(z, V - DV_STEP, T, Pin, dPdV=dPdV_minus)

      f2 = (dPdV_plus - dPdV_minus) / (2.0_pr * DV_STEP)
   end subroutine residual_critical_eqs


   ! ========================================================================
   ! Jacobian: matrix of derivatives ∂f_i / ∂x_j
   ! ========================================================================
   subroutine jacobian_critical_eqs(model, nc, i_comp, V, T, J)
      !! # jacobian_critical_eqs
      !! Approximate the Jacobian of \((P_V,P_{VV})\) with respect to \((V,T)\).
      !!
      !! # Description
      !! Centered finite differences use the module increments `DV_STEP` and
      !! `DT_STEP`.
      class(ArModel), intent(in) :: model !! Equation of state used to evaluate the critical residual Jacobian.
      integer, intent(in) :: i_comp !! One-based index of the pure component.
      integer, intent(in) :: nc !! Number of components.
      real(pr), intent(in) :: V !! Total volume [L].
      real(pr), intent(in) :: T !! Temperature [K].
      real(pr), intent(out) :: J(2,2) !! Jacobian rows \((P_V,P_{VV})\), columns \((V,T)\).

      real(pr) :: f1_V_plus
      real(pr) :: f1_V_minus
      real(pr) :: f1_T_plus
      real(pr) :: f1_T_minus
      real(pr) :: f2_V_plus
      real(pr) :: f2_V_minus
      real(pr) :: f2_T_plus
      real(pr) :: f2_T_minus
      real(pr) :: P_dummy

      ! ∂f1/∂V and ∂f2/∂V
      call residual_critical_eqs(model, nc, i_comp, V + DV_STEP, T, f1_V_plus, f2_V_plus, P_dummy)
      call residual_critical_eqs(model, nc, i_comp, V - DV_STEP, T, f1_V_minus, f2_V_minus, P_dummy)
      J(1,1) = (f1_V_plus - f1_V_minus) / (2.0_pr * DV_STEP)
      J(2,1) = (f2_V_plus - f2_V_minus) / (2.0_pr * DV_STEP)

      ! ∂f1/∂T and ∂f2/∂T
      call residual_critical_eqs(model, nc, i_comp, V, T + DT_STEP, f1_T_plus, f2_T_plus, P_dummy)
      call residual_critical_eqs(model, nc, i_comp, V, T - DT_STEP, f1_T_minus, f2_T_minus, P_dummy)
      J(1,2) = (f1_T_plus - f1_T_minus) / (2.0_pr * DT_STEP)
      J(2,2) = (f2_T_plus - f2_T_minus) / (2.0_pr * DT_STEP)
   end subroutine jacobian_critical_eqs

end module yaeos__critical_pure_point_solver
