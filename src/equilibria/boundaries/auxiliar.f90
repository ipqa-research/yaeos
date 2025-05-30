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

   subroutine detect_critical(nc, np, kinds_x, kind_w, binary_stop, X, dXdS, ns, dS, S)
      !! # detect_critical
      !! Detect if the system is close to a critical point.
      !!
      !! # Description
      !! When the system is close to a critical point, the \(\ln K_i^l\) values
      !! are close to zero, since the composition of the incipient phase and the
      !! \(l\) phase are similar (equal in the critical point). This can be used
      !! to detect if the system is close to a critical point and force a jump
      !! above it.
      !!
      !! # References
      !!
      integer, intent(in) :: nc
      !! Number of components in the mixture.
      integer, intent(in) :: np
      !! Number of main phases.
      character(len=14), intent(in out) :: kinds_x(np)
      !! Kinds of the main phases.
      character(len=14), intent(in out) :: kind_w
      !! Kind of the incipient phase.
      logical, intent(in ) :: binary_stop
      !! If true, stop at the critical point if its a binary system.
      real(pr), intent(in out) :: X(:)
      !! Vector of variables.
      real(pr), intent(in out) :: dXdS(:)
      !! Sensitivity of the variables wrt the specification.
      integer, intent(in out) :: ns
      !! Number of the specified variable.
      real(pr), intent(in out) :: dS
      !! Step size of the specification for the next point.
      real(pr), intent(in out) :: S
      !! Specification value.

      real(pr) :: Xold(size(X))
      character(len=14) :: incipient_kind
      integer :: i, lb, ub

      Xold = X

      do i=1,np
         lb = (i-1)*nc + 1
         ub = i*nc

         do while(maxval(abs(X(lb:ub))) < 0.01)
            if (nc == 2 .and. maxval(abs(X(lb:ub))) < 1e-6) then
               dS=0
               return
            end if
            X = X + dXdS * dS
         end do

         if (all(Xold(lb:ub) * (X(lb:ub) + dXdS(lb:ub)*dS) < 0)) then
            incipient_kind = kind_w
            kind_w = kinds_x(i)
            kinds_x(i) = incipient_kind
            ! Interpolate here

            if (nc == 2 .and. binary_stop) then
               dS=0
               return
            end if
         end if

      end do

   end subroutine detect_critical
end module yaeos__equilibria_boundaries_auxiliar
