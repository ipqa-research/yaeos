module yaeos__critical_pure_point_solver
   !! # `critical_pure_point_solver`
   !! Calculation of pure-component critical points
   !! using the yaeos equation of state and numerical derivatives.
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

   use yaeos__constants, only: pr
   use yaeos__models_ar, only: ArModel

   implicit none
   private
   public :: find_critical_points_all_components

   integer, parameter :: MAX_ITER = 151
   real(pr), parameter :: TOL_RES = 2.0e-8_pr      ! residual tolerance
   real(pr), parameter :: TOL_STEP = 2.0e-6_pr     ! step tolerance

   ! Steps for numerical derivatives
   real(pr), parameter :: DV_STEP = 2.0e-4_pr      ! L/mol
   real(pr), parameter :: DT_STEP = 2.0e-2_pr      ! K

contains

   ! ========================================================================
   ! MAIN ROUTINE: Find critical point for each pure component
   ! ========================================================================
   subroutine find_critical_points_all_components(model, nc, Vc, Tc, Pc, converged)
      class(ArModel), intent(in) :: model
      integer, intent(in) :: nc
      real(pr), intent(out) :: Vc(nc)
      real(pr), intent(out) :: Tc(nc)
      real(pr), intent(out) :: Pc(nc)
      logical, intent(out) :: converged(nc)

      integer :: i
      real(pr) :: V_init, T_init
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
      class(ArModel), intent(in) :: model
      integer, intent(in) :: i
      integer, intent(in) :: nc
      real(pr), intent(out) :: V_est
      real(pr), intent(out) :: T_est

      real(pr) :: z(nc), T, P_sat, VL, VV, ratio
      logical :: psat_converged
      integer :: iter

      z = 0.0_pr
      z(i) = 1.0_pr

      ! Initial temperature heuristic scan
      T = 50.0_pr
      V_est = 0.15_pr
      T_est = T

      do iter = 1, 20
         P_sat = model%Psat_pure(i, T, Vl=VL, Vv=Vv, converged=psat_converged)
         if (VL > 0.0_pr .and. VV > VL) then
            ratio = VV / VL
            ! Cailletet-Mathias mean rectilinear volume estimate
            V_est = (2.0_pr * VL * VV) / (VL + VV)

            if (ratio < 1.25_pr) then
               T_est = T
               return
            end if

            ! Advance temperature upward towards Tc
            T = T * (1.0_pr + 0.08_pr * log(ratio))
         else
            ! Step temperature down if above pseudo-critical or single phase
            T = T * 0.75_pr
         end if
      end do
      T_est = max(T, 50.0_pr)
   end subroutine estimate_critical_point


   ! ========================================================================
   ! SOLVER: Newton-Raphson 2D to solve [dPdV, d²PdV²] = 0
   ! ========================================================================
   subroutine solve_critical_point(model, nc, i_comp, V_init, T_init, Vc, Tc, Pc, converged)
      class(ArModel), intent(in) :: model
      integer, intent(in) :: nc, i_comp
      real(pr), intent(in) :: V_init, T_init
      real(pr), intent(out) :: Vc, Tc, Pc
      logical, intent(out) :: converged

      real(pr) :: V, T, P
      real(pr) :: f1, f2, norm_f
      real(pr) :: J(2,2), J_inv(2,2), det, dX(2)
      integer :: iter

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
      class(ArModel), intent(in) :: model
      integer, intent(in) :: i_comp, nc
      real(pr), intent(in) :: V, T
      real(pr), intent(out) :: f1, f2, P

      real(pr) :: z(nc), dPdV, dPdV_plus, dPdV_minus, Pin

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
      class(ArModel), intent(in) :: model
      integer, intent(in) :: i_comp, nc
      real(pr), intent(in) :: V, T
      real(pr), intent(out) :: J(2,2)

      real(pr) :: f1_V_plus, f1_V_minus, f1_T_plus, f1_T_minus
      real(pr) :: f2_V_plus, f2_V_minus, f2_T_plus, f2_T_minus
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
