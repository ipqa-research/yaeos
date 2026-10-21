module yaeos__equilibria_flash_second_order
   !! # Second-order isothermal flash
   !!
   !! ## Description
   !! Two-phase flash at specified temperature and pressure by direct
   !! minimization of the reduced Gibbs energy, following the formulation of
   !! Michelsen (1982) described in `AGENTS.md`.
   !!
   !! The minimization variables are the scaled heavy-phase mole numbers
   !! \(\theta_i = l_i/z_i\), so the light phase is \(v_i = z_i(1-\theta_i)\)
   !! and the admissible region is \(0 < \theta_i < 1\). The Gibbs energy is
   !! minimized with a Newton step based on a modified Cholesky decomposition
   !! (Murray), combined with a line search that first computes the largest
   !! step compatible with the bounds and then backtracks on the Gibbs energy.
   !!
   !! For `ArModel` the heavy phase is always evaluated with the liquid root.
   !! The light phase uses the stable root by default, but it can be forced to
   !! the liquid root with `light_root="liquid"` (useful for liquid-liquid
   !! splits). For `GeModel`, both phases are described by activity
   !! coefficients and the initial split is obtained from stability analysis.
   use yaeos__constants, only: pr
   use yaeos__models_base, only: BaseModel
   use yaeos__models_ar, only: ArModel
   use yaeos__models_ge, only: GeModel
   use yaeos__equilibria_equilibrium_state, only: EquilibriumState
   use yaeos__equilibria_rachford_rice, only: betato01
   use yaeos__equilibria_auxiliar, only: k_wilson
   use yaeos__equilibria_stability, only: min_tpd
   implicit none

   private
   public :: flash_second_order

contains

   type(EquilibriumState) function flash_second_order( &
      model, z, T, P_spec, k0, light_root, iters, tol, max_iters &
      ) result(flash)
      !! # `flash_second_order`
      !!
      !! ## Description
      !! Find a two-phase split at fixed temperature and pressure by minimizing
      !! the reduced Gibbs energy with a second-order method.
      !!
      !! ## Arguments
      !! - `model`: thermodynamic model (`ArModel` or `GeModel`).
      !! - `z`: global mole fractions.
      !! - `T`: temperature [K].
      !! - `P_spec`: pressure [bar].
      !! - `k0`: optional initial K-factors. If omitted, Wilson factors are
      !!   used for `ArModel` and a stability analysis is used for `GeModel`.
      !! - `light_root`: root used for the light phase with `ArModel`
      !!   (`"stable"` by default, `"liquid"` for a liquid-liquid split).
      !! - `iters`: number of iterations used.
      !! - `tol`: gradient-norm tolerance (default `1e-8`).
      !! - `max_iters`: maximum number of iterations (default `200`).
      class(BaseModel), intent(in) :: model !! Thermodynamic model
      real(pr), intent(in) :: z(:) !! Global mole fractions
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: P_spec !! Pressure [bar]
      real(pr), optional, intent(in) :: k0(:) !! Initial K-factors (y/x)
      character(len=*), optional, intent(in) :: light_root
      !! Root of the light phase for `ArModel` (`"stable"` or `"liquid"`)
      integer, optional, intent(out) :: iters !! Iterations used
      real(pr), optional, intent(in) :: tol !! Convergence tolerance
      integer, optional, intent(in) :: max_iters !! Maximum iterations

      real(pr) :: tolerance, grad_norm
      integer :: max_iterations, i, iterations

      real(pr) :: theta(size(z)), theta_trial(size(z))
      real(pr) :: l(size(z)), v(size(z)), x(size(z)), y(size(z))
      real(pr) :: K(size(z)), w(size(z))
      real(pr) :: Ltot, Vtot, beta
      real(pr) :: mintpd

      real(pr) :: Gibbs, g_l(size(z)), H_l(size(z), size(z))
      real(pr) :: g(size(z)), H(size(z), size(z))
      real(pr) :: step(size(z))
      real(pr) :: alpha, alpha_max, G_trial, tolG
      logical :: accepted
      integer :: j

      character(len=:), allocatable :: light_root_type

      tolerance = 1e-8_pr
      if (present(tol)) tolerance = tol
      max_iterations = 200
      if (present(max_iters)) max_iterations = max_iters
      light_root_type = "stable"
      if (present(light_root)) light_root_type = trim(light_root)

      ! ========================================================================
      ! Initialization
      ! ------------------------------------------------------------------------
      if (present(k0)) then
         K = k0
      else
         select type (model)
          class is (ArModel)
            K = k_wilson(model, T, P_spec)
          class is (GeModel)
            ! A stability analysis provides the composition of the incipient
            ! phase. It is the natural initial estimate when no vapour-like
            ! root exists in the model.
            call min_tpd(model, z, P_spec, T, mintpd, w)
            K = w/z
         end select
      end if

      ! A null K-factor would put all of a component in a single phase and make
      ! the composition logarithms undefined. A tiny positive value keeps the
      ! initialization inside the admissible region.
      K = max(K, 1e-30_pr)

      ! Make sure that the K-factors are usable for a two-phase initialization.
      call betato01(z, K)

      ! Start from a balanced split. The equilibrium phase fraction is usually
      ! far from zero, and using a balanced initial amount avoids starting the
      ! minimization with one of the phases almost absent.
      beta = 0.5_pr
      l = z*(1 - beta)/(1 + beta*(K - 1))
      theta = l/z

      ! ========================================================================
      ! Gibbs energy minimization
      ! ------------------------------------------------------------------------
      ! The minimization is carried out in the scaled variables
      ! \f$\theta_i = l_i/z_i\f$, so the admissible region is simply
      ! \f$0 < \theta_i < 1\f$.
      iterations = 0
      do while (iterations < max_iterations)
         l = z*theta
         call evaluate_gibbs( &
            model, z, l, T, P_spec, light_root_type, Gibbs, g_l, H_l &
            )

         ! Chain rule: g_theta = Z g_l, H_theta = Z H_l Z.
         g = z*g_l
         do j = 1, size(z)
            H(:, j) = z*H_l(:, j)*z(j)
         end do

         grad_norm = maxval(abs(g))
         if (grad_norm < tolerance) exit

         iterations = iterations + 1

         ! Modified Cholesky: H + E is positive definite, so `step` is a
         ! descent direction for the Gibbs energy.
         call solve_modified_cholesky(H, -g, step)

         ! Largest step that keeps 0 < theta < 1 in every component.
         alpha_max = huge(1.0_pr)
         do j = 1, size(z)
            if (step(j) > 0.0_pr) then
               alpha_max = min(alpha_max, (1.0_pr - theta(j))/step(j))
            else if (step(j) < 0.0_pr) then
               alpha_max = min(alpha_max, -theta(j)/step(j))
            end if
         end do
         if (alpha_max < huge(1.0_pr)) alpha_max = 0.99_pr*alpha_max

         ! Backtracking line search on the Gibbs energy. The acceptance allows
         ! an increase of `tolG`, which absorbs the evaluation noise of the EoS
         ! volume solver (~1e-12 for cubic models). Without this margin the
         ! line search rejects the full Newton step at the noise floor and the
         ! iteration stalls while the gradient is still around 1e-8.
         tolG = 1.0e4_pr*epsilon(1.0_pr)*max(1.0_pr, abs(Gibbs))
         alpha = min(1.0_pr, alpha_max)
         G_trial = Gibbs
         accepted = .false.
         do i = 1, 60
            theta_trial = theta + alpha*step
            if (all(theta_trial > 0 .and. theta_trial < 1)) then
               l = z*theta_trial
               call evaluate_gibbs( &
                  model, z, l, T, P_spec, light_root_type, &
                  G_trial, g_l, H_l &
                  )
               if (G_trial <= Gibbs + tolG) then
                  accepted = .true.
                  exit
               end if
            end if
            alpha = 0.5_pr*alpha
            if (alpha < 1e-16_pr) exit
         end do

         if (.not. accepted) exit
         theta = theta_trial
      end do

      ! ========================================================================
      ! Format results
      ! ------------------------------------------------------------------------
      ! Re-evaluate at the final point so the convergence test always uses the
      ! gradient of the phase split that is actually returned.
      l = z*theta
      call evaluate_gibbs( &
         model, z, l, T, P_spec, light_root_type, Gibbs, g_l, H_l &
         )
      g = z*g_l

      v = z - l
      Ltot = sum(l)
      Vtot = sum(v)
      x = l/Ltot
      y = v/Vtot

      flash%kind = "split"
      flash%T = T
      flash%P = P_spec
      flash%x = x
      flash%y = y
      flash%beta = Vtot
      flash%iters = iterations

      ! A small gradient and two different phases are required. This rejects
      ! non-converged calculations and the trivial solution x = y.
      if (maxval(abs(g)) > tolerance .or. maxval(abs(x - y)) < 1e-10_pr) then
         flash%kind = "failed"
         flash%P = -1.0_pr
         if (present(iters)) iters = iterations
         return
      end if

      select type (model)
       class is (ArModel)
         call model%volume(x, P_spec, T, flash%Vx, root_type="liquid")
         call model%volume(y, P_spec, T, flash%Vy, root_type=light_root_type)
       class is (GeModel)
         flash%Vx = 0.0_pr
         flash%Vy = 0.0_pr
      end select

      if (present(iters)) iters = iterations
   end function flash_second_order

   subroutine evaluate_gibbs(model, z, l, T, P, light_root, Gibbs, g, H)
      !! Evaluate the reduced Gibbs energy, its gradient and its Hessian with
      !! respect to the heavy-phase mole numbers.
      class(BaseModel), intent(in) :: model
      real(pr), intent(in) :: z(:) !! Global mole numbers (1 mol feed)
      real(pr), intent(in) :: l(:) !! Heavy-phase mole numbers
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: P !! Pressure [bar]
      character(len=*), intent(in) :: light_root
      real(pr), intent(out) :: Gibbs !! \(G/RT\)
      real(pr), intent(out) :: g(:) !! Gradient
      real(pr), intent(out) :: H(:, :) !! Hessian

      integer :: i, j, nc
      real(pr) :: v(size(z)), x(size(z)), y(size(z))
      real(pr) :: Ltot, Vtot
      real(pr) :: lnphi_l(size(z)), lnphi_v(size(z))
      real(pr) :: dlnphi_l(size(z), size(z)), dlnphi_v(size(z), size(z))

      nc = size(z)
      v = z - l
      Ltot = sum(l)
      Vtot = sum(v)
      x = l/Ltot
      y = v/Vtot

      select type (model)
       class is (ArModel)
         call model%lnphi_pt( &
            l, P, T, root_type="liquid", &
            lnPhi=lnphi_l, dlnPhidn=dlnphi_l &
            )
         call model%lnphi_pt( &
            v, P, T, root_type=light_root, &
            lnPhi=lnphi_v, dlnPhidn=dlnphi_v &
            )
       class is (GeModel)
         call model%ln_activity_coefficient( &
            l, T, lngamma=lnphi_l, dlngammadn=dlnphi_l &
            )
         call model%ln_activity_coefficient( &
            v, T, lngamma=lnphi_v, dlngammadn=dlnphi_v &
            )
      end select

      Gibbs = sum(l*(log(x) + lnphi_l)) + sum(v*(log(y) + lnphi_v))
      g = log(x) + lnphi_l - log(y) - lnphi_v

      ! Q + Phi, with Q the ideal-mixing contribution and Phi the fugacity
      ! coefficient derivatives of both phases.
      do j = 1, nc
         do i = 1, nc
            H(i, j) = -1/Ltot - 1/Vtot + dlnphi_l(i, j) + dlnphi_v(i, j)
         end do
         H(j, j) = H(j, j) + 1/l(j) + 1/v(j)
      end do
   end subroutine evaluate_gibbs

   subroutine solve_modified_cholesky(H, b, x)
      !! Solve (H + D) x = b with Murray's modified Cholesky (1972), Sec. 4.9.
      !! Factorise by columns H + D = Lb*Lb^T with |l(j,i)| <= beta (Eqs. 4.9.7-12).
      !! If H is sufficiently positive definite then D = 0 (pure Newton).
      real(pr), intent(in)  :: H(:, :)
      real(pr), intent(in)  :: b(:)
      real(pr), intent(out) :: x(:)

      integer  :: i, j, n
      real(pr) :: gamma, xi, beta, delta, theta, lh
      real(pr) :: Lb(size(H, 1), size(H, 1)), lj(size(H, 1)), y(size(H, 1))

      n = size(H, 1)

      ! Ec. 4.9.16
      gamma = 0.0_pr
      xi = 0.0_pr
      do i = 1, n
         gamma = max(gamma, abs(H(i, i)))
         do j = 1, n
            if (j /= i) xi = max(xi, abs(H(i, j)))
         end do
      end do
      beta = sqrt(max(gamma, xi/max(1.0_pr, sqrt(real(n, pr)**2 - 1.0_pr)), &
         epsilon(1.0_pr)))

      delta = sqrt(epsilon(1.0_pr))

      Lb = 0.0_pr
      do i = 1, n
         ! Ec. 4.9.8
         lh = max(delta, sqrt(abs(H(i, i) - sum(Lb(i, 1:i - 1)**2))))

         ! Ec. 4.9.9
         theta = 0.0_pr
         do j = i + 1, n
            lj(j) = (H(j, i) - dot_product(Lb(i, 1:i - 1), Lb(j, 1:i - 1)))/lh
            theta = max(theta, abs(lj(j)))
         end do

         if (theta <= beta) then
            ! Ec. 4.9.10
            Lb(i, i) = lh
            do j = i + 1, n
               Lb(j, i) = lj(j)
            end do
         else
            ! Ec. 4.9.11-12
            Lb(i, i) = theta*lh/beta
            do j = i + 1, n
               Lb(j, i) = beta*lj(j)/theta
            end do
         end if
      end do

      ! Lb y = b
      do i = 1, n
         y(i) = (b(i) - dot_product(Lb(i, 1:i - 1), y(1:i - 1)))/Lb(i, i)
      end do
      ! Lb^T x = y
      do i = n, 1, -1
         x(i) = (y(i) - dot_product(Lb(i + 1:n, i), x(i + 1:n)))/Lb(i, i)
      end do
   end subroutine solve_modified_cholesky
end module yaeos__equilibria_flash_second_order
