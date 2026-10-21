module yaeos__equilibria_flash_second_order
   !! # Second-order isothermal flash
   !!
   !! ## Description
   !! Two-phase flash at specified temperature and pressure by direct
   !! minimization of the reduced Gibbs energy, following the formulation of
   !! Michelsen (1982) described in `AGENTS.md`.
   !!
   !! The independent variables are the mole numbers of the heavy phase,
   !! \(l_i\). The light-phase mole numbers follow from the material balance
   !! \(v_i = z_i - l_i\). The Gibbs energy is minimized with a Newton step
   !! based on a modified Cholesky decomposition (Murray), combined with a
   !! backtracking line search that keeps the iterate inside the physically
   !! admissible region \(0 < l_i < z_i\).
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

      real(pr) :: l(size(z)), v(size(z)), x(size(z)), y(size(z))
      real(pr) :: K(size(z)), w(size(z))
      real(pr) :: Ltot, Vtot, beta
      real(pr) :: mintpd

      real(pr) :: Gibbs, g(size(z)), H(size(z), size(z))
      real(pr) :: step(size(z)), l_trial(size(z))
      real(pr) :: alpha, G_trial
      logical :: accepted

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

      ! ========================================================================
      ! Gibbs energy minimization
      ! ------------------------------------------------------------------------
      iterations = 0
      do while (iterations < max_iterations)
         call evaluate_gibbs( &
            model, z, l, T, P_spec, light_root_type, Gibbs, g, H &
            )

         grad_norm = maxval(abs(g))
         if (grad_norm < tolerance) exit

         iterations = iterations + 1

         ! Modified Cholesky: H + E is positive definite, so `step` is a
         ! descent direction for the Gibbs energy.
         call solve_modified_cholesky(H, -g, step)

         ! A trust-region cap avoids extremely large steps when the Hessian
         ! is ill-conditioned near a phase boundary.
         if (maxval(abs(step)) > 0.5_pr*maxval(z)) then
            step = step*(0.5_pr*maxval(z)/maxval(abs(step)))
         end if

         ! Backtracking line search. The physical bounds also prevent the
         ! iterate from reaching the trivial solution (l = 0 or v = 0).
         alpha = 1.0_pr
         G_trial = Gibbs
         accepted = .false.
         do i = 1, 40
            l_trial = l + alpha*step
            if (all(l_trial > 0 .and. l_trial < z)) then
               call evaluate_gibbs( &
                  model, z, l_trial, T, P_spec, light_root_type, &
                  G_trial, g, H &
                  )
               if (G_trial <= Gibbs) then
                  accepted = .true.
                  exit
               end if
            end if
            alpha = 0.5_pr*alpha
         end do

         if (.not. accepted) exit
         l = l_trial
      end do

      ! ========================================================================
      ! Format results
      ! ------------------------------------------------------------------------
      ! Re-evaluate at the final point so the convergence test always uses the
      ! gradient of the phase split that is actually returned.
      call evaluate_gibbs( &
         model, z, l, T, P_spec, light_root_type, Gibbs, g, H &
         )

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
      !! Solve \((H + E)x = b\) using a modified Cholesky factorisation.
      !!
      !! The diagonal pivots are bounded from below. This is equivalent to
      !! adding a diagonal matrix \(E\) to \(H\), and guarantees that the
      !! resulting matrix is positive definite. If \(H\) is already positive
      !! definite (with comfortably positive pivots) then \(E = 0\).
      real(pr), intent(in) :: H(:, :)
      real(pr), intent(in) :: b(:)
      real(pr), intent(out) :: x(:)

      integer :: i, j, k, n
      real(pr) :: pivot, min_pivot, scale, value
      real(pr) :: Lmat(size(H, 1), size(H, 2))
      real(pr) :: D(size(H, 1)), y(size(H, 1))

      n = size(H, 1)
      scale = max(1.0_pr, maxval(abs(H)))
      min_pivot = 1e-12_pr*scale

      Lmat = 0.0_pr
      do i = 1, n
         Lmat(i, i) = 1.0_pr
      end do
      D = 0.0_pr

      do j = 1, n
         pivot = H(j, j)
         do k = 1, j - 1
            pivot = pivot - Lmat(j, k)**2*D(k)
         end do
         if (pivot < min_pivot) pivot = min_pivot
         D(j) = pivot

         do i = j + 1, n
            value = H(i, j)
            do k = 1, j - 1
               value = value - Lmat(i, k)*Lmat(j, k)*D(k)
            end do
            Lmat(i, j) = value/pivot
         end do
      end do

      ! Solve L y = b, D z = y, L^T x = z.
      do i = 1, n
         y(i) = b(i)
         do k = 1, i - 1
            y(i) = y(i) - Lmat(i, k)*y(k)
         end do
      end do
      do i = 1, n
         y(i) = y(i)/D(i)
      end do
      do i = n, 1, -1
         x(i) = y(i)
         do k = i + 1, n
            x(i) = x(i) - Lmat(k, i)*x(k)
         end do
      end do
   end subroutine solve_modified_cholesky
end module yaeos__equilibria_flash_second_order
