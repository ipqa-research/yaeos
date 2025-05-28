module yaeos__equilibria_boundaries_auxiliar
   !! Equilibria boundaries auxiliar module
   !! This module contains the auxiliar functions and subroutines
   !! used in the phase-boundaries calculations.
   use yaeos__constants, only: R, pr
   implicit none
contains
   subroutine get_z(alpha, z_0, z_inj, z, dzda)
      !! Calculate the fluid composition based on an amount of addition
      !! of second fluid.
      !!
      !! The injection can be considered as two kinds of injection:
      !! - Displacement: \( z = \alpha z_i + (1-\alpha) z_0 \)
      !! - Addition:  \( z = \frac{\alpha z_i + (1-\alpha) z_0}{\sum_{i=1}^N \alpha z_i + (1-\alpha) z_0} \)
      real(pr), intent(in)  :: alpha !! Addition percentaje \( \alpha \)
      real(pr), intent(in) :: z_inj(:)
      real(pr), intent(in) :: z_0(:)
      real(pr), intent(out) :: z(size(z_0)) !! New composition
      real(pr), optional, intent(out) :: dzda(size(z_0)) !! Derivative wrt \(\alpha\)

      z = z_inj * alpha + (1.0_pr - alpha)*z_0
      if (present(dzda)) dzda = z_inj - z_0
   end subroutine get_z
end module yaeos__equilibria_boundaries_auxiliar
