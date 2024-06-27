module yaeos__phase_equilibria_rachford_rice
   use yaeos__constants, only: pr
   implicit none
contains
   subroutine betato01(z, K)
      !! Modify K-factor values to assure that \(\beta\) lies between (0,1)
      implicit none
      real(pr), intent(in) :: z(:) !! Molar fractions of the system
      real(pr) :: K(:) !! K factors \(\frac{y_i}{x_i}\)

      real(pr) :: g0, g1  ! function g valuated at beta=0 and 1, based on K factors

      g1 = 1.0
      do while (g0 < 0 .or. g1 > 0)
         g0 = sum(z*K) - 1._pr
         g1 = 1._pr - sum(z/K)
         if (g0 < 0) then
            ! Increased volatiliy will bring the solution from
            ! subcooled liquid into VLE
            K = 1.1_pr * K
         else if (g1 > 0) then
            ! Decreased volatiliy will bring the solution from
            ! superheated vapor into VLE
            K = 0.9_pr * K
         end if
      end do
   end subroutine betato01

   subroutine betalimits(z, K, bmin, bmax)
      !! Define beta limits to avoid overshooting when solving the Rachford-Rice
      !! equation.
      !!
      !! This is based on the assumtion that either \(y_i < 1\) and \(x_i < 1\).
      real(pr), intent(in) :: z(:) !! Molar fractions vector
      real(pr), intent(in) :: K(:) !! K-factors
      real(pr), intent(out) :: bmin !! Minimum beta value
      real(pr), intent(out) :: bmax  !! Maximum beta value

      real(pr), dimension(size(z)) :: vmin, vmax

      vmin = 0.d0
      ! max=1.001d0    ! modified  3/3/15 (not to generate false separations with beta 0.9999...)
      vmax = 1.00001_pr ! modified 28/6/15 (to prevent overshooting in the Newton for solving RR eq.)

      where (K * z > 1)
         vmin = (K * z - 1._pr) / (K - 1._pr)
      elsewhere (K < z)
         vmax = (1 - z)/(1 - K)
      end where

      bmin = maxval(vmin)
      bmax = minval(vmax)
   end subroutine betalimits

   subroutine rachford_rice(z, K, beta, rr, drrdb)
      !! Rachford-Rice equation for a two phase system
      real(pr), intent(in) :: z(:)
      real(pr), intent(in) :: K(:)
      real(pr), intent(in) :: beta

      real(pr), intent(out) :: rr
      real(pr), intent(out) :: drrdb

      real(pr) :: denom(size(z))

      denom = 1 + beta*(K - 1._pr)
      rr = sum(z*(K - 1._pr)/denom)
      drrdb = -sum(z*(K - 1._pr)**2/denom**2)
   end subroutine rachford_rice

   subroutine solve_rr(z, K, beta, beta_min, beta_max)
      !! Solve the Rachford-Rice Equation.
      real(pr), intent(in) :: z(:) !! Mole fractions vector
      real(pr), intent(in) :: K(:) !! K-factors
      real(pr), intent(out) :: beta_min !! Lower limit for \(\beta\)
      real(pr), intent(out) :: beta_max !! Upper limit for \(\beta\)
      real(pr), intent(out) :: beta !! \(\beta\) value

      real(pr) :: g, dgdb
      real(pr) :: step

      g = 1.0
      step = 1.0

      call betalimits(z, k, beta_min, beta_max)

      do while (abs(g) > 1.d-5 .and. abs(step) > 1.d-10)
         call rachford_rice(z, k, beta, g, dgdb)
         step = -g/dgdb
         beta = beta + step
         do while ((beta < beta_min .or. beta_max < beta) .and. step > 1e-10)
            step = step/2
            beta = beta - step
         end do
      end do
   end subroutine solve_rr
end module yaeos__phase_equilibria_rachford_rice
