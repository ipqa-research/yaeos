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

   subroutine detect_critical(&
         nc, np, point, kinds_x, kind_w, binary_stop, Xold, X, dXdS, ns, dS, S, found_critical, Xc&
      )
      !! # detect_critical
      !! Detect if the system is close to a critical point.
      !!
      !! # Description
      !! When the system is close to a critical point, the \(\ln K_i^l\) values
      !! are close to zero, since the composition of the incipient phase and the
      !! \(l\) phase are similar (equal in the critical point). This can be used
      !! to detect if the system is close to a critical point and force a jump
      !! above it.
      integer, intent(in) :: nc
      !! Number of components in the mixture.
      integer, intent(in) :: np
      !! Number of main phases.
      integer, intent(in) :: point
      !! Point number in the phase boundary.
      character(len=14), intent(in out) :: kinds_x(np)
      !! Kinds of the main phases.
      character(len=14), intent(in out) :: kind_w
      !! Kind of the incipient phase.
      logical, intent(in) :: binary_stop
      !! If true, stop at the critical point if its a binary system.
      real(pr), intent(in) :: Xold(:)
      !! Old vector of variables.
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
      logical, intent(out) :: found_critical
      !! If true, a critical point was found.
      real(pr) :: Xc(size(X))
      !! Vector of variables at the critical point.
      
      character(len=14) :: incipient_kind
      integer :: i, lb, ub
      integer :: ncomp
      real(pr) :: a,  Xnew(size(X))

      real(pr) :: lnKold(nc), lnK(nc)
         
      found_critical = .false.

      do i=1,np
         lb = (i-1)*nc + 1
         ub = i*nc
         ! TODO: In many cases this makes more damage than good.
         ! do while(maxval(abs(X(lb:ub))) < 0.05)
         !    print *, "Increasing X", point, X(lb:ub)
         !    dS = sign(max(abs(dS), sqrt(abs(X(ns))/10), 0.01), dS)
         !    if (nc == 2 .and. maxval(abs(X(lb:ub))) < 1e-6 .and. binary_stop) then
         !       ! Reached to a critical point in a Txy/Pxy calculation for a
         !       ! binary system, stop the calculation.
         !       dS=0
         !       return
         !    end if
         !    X = X + dXdS * dS
         !    print *, "next: X", point, X(lb:ub)
         ! end do

         lnKold = Xold(lb:ub)
         lnK = X(lb:ub) + dXdS(lb:ub) * dS
         Xnew = X + dXdS * dS

         if (point > 1 .and. all(lnKold * lnK < 0)) then
            ! In Liquid-Liquid lines that start from a critical point, this
            ! could be a false positive, so we check that the point is not
            ! the first one.
            found_critical = .true.
            incipient_kind = kind_w
            kind_w = kinds_x(i)
            kinds_x(i) = incipient_kind
            
            ! 0 = a*Xnew(ns) + (1-a)*X(ns) < Interpolation equation to get X(ns) = 0
            ncomp = maxloc(abs(lnK - lnKold), dim=1)
            a = -lnKold(ncomp)/(lnK(ncomp) - lnKold(ncomp))
            Xc = a * Xnew + (1-a)*Xold

            if (nc == 2 .and. binary_stop) then
               dS=0
               return
            end if

            dS = sign(max(abs(dS), sqrt(abs(X(ns))/10), 0.01), dS)
            X = Xc
            return
         end if
      end do
   end subroutine detect_critical
end module yaeos__equilibria_boundaries_auxiliar
