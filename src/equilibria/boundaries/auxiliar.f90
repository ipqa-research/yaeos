module yaeos__equilibria_boundaries_auxiliar
   !! Equilibria boundaries auxiliar module
   !! This module contains the auxiliar functions and subroutines
   !! used in the phase-boundaries calculations.
   use yaeos__constants, only: R, pr
   use yaeos__math, only: interpol
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
      nc, np, point, kinds_x, kind_w, binary_stop, Xold, &
      X, dXdS, ns, dS, S, found_critical, Xc&
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
      use yaeos__math, only: interpol
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

      real(pr) :: T, P

      real(pr) :: limit

      T = exp(X(np*nc+np+2))
      P = exp(X(np*nc+np+1))

      found_critical = .false.

      limit = 0.01 + nc * (0.1 - 0.01)/(20. - 2)

      do i=1,np
         lb = (i-1)*nc + 1
         ub = i*nc
         ! TODO: In many cases this makes more damage than good.
         do while(maxval(abs(X(lb:ub))) < min(limit, 0.05_pr))
            if (nc == 2 .and. maxval(abs(X(lb:ub))) < 1e-6 .and. binary_stop) then
               ! Reached to a critical point in a Txy/Pxy calculation for a
               ! binary system, stop the calculation.
               dS=0
               return
            end if

            if (maxval(abs(dXdS(lb:ub)*dS)) > 0.07) then
               X = X - dXdS * dS / 2._pr
            else
               X = X + dXdS * dS
            end if
         end do

         Xnew = X + dXdS * dS

         lnKold = Xold(lb:ub)
         lnK = Xnew(lb:ub)

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
            Xc = interpol(lnKold(ncomp), lnK(ncomp), Xold, Xnew, 0.0_pr)

            if (nc == 2 .and. binary_stop) then
               dS=0
               return
            end if

            ! Start from the critical point and then do small steps until
            ! we are a bit far from it.
            X = Xc + sign(0.001_pr, dS) * dXdS
            ! do while(near_critical(nc, np, X))
            do while(maxval(abs(X(lb:ub))) < maxval(abs(Xold(lb:ub))))
               X = X + sign(0.001_pr, dS) * dXdS
            end do

            if (ns > 0) then
               S = X(ns)
            end if
            return
         end if
      end do
   end subroutine detect_critical

   subroutine check_critical_jump(nc, np, ns, kinds_x, kind_w, X, X_last_converged, Xc, found_critical, jumped_critical)
      use yaeos__math, only: interpol
      integer, intent(in) :: nc !! Number of components
      integer, intent(in) :: np !! Number of main phases
      integer, intent(in) :: ns !! Number of the specified variable
      character(len=14), intent(inout) :: kinds_x(np) !! Kinds of the main phases
      character(len=14), intent(inout) :: kind_w !! Kind of the incipient phase
      real(pr), intent(in) :: X(:) !! Current point
      real(pr), intent(in) :: X_last_converged(:) !! Previously converged point
      real(pr), intent(out) :: Xc(:) !! Critical point
      logical, intent(in) :: found_critical !! If a critical point was already found in the current point
      logical, intent(out) :: jumped_critical !! If a critical point was jumped

      character(len=14) :: incipient_kind

      integer :: l, lb, ub
      integer :: ncomp
      jumped_critical = .false.
      do l=1,np
         lb = (l-1)*nc + 1
         ub = l*nc
         if (all(X(lb:ub) * X_last_converged(lb:ub) < 0._pr)) then
            ncomp = maxloc(abs(X(lb:ub) - X_last_converged(lb:ub)), dim=1) + lb - 1
            jumped_critical = .true.
            Xc = interpol(X_last_converged(ncomp), X(ncomp), X_last_converged, X, 0._pr)

            if (jumped_critical .and. .not. found_critical) then
               incipient_kind = kind_w
               kind_w = kinds_x(l)
               kinds_x(l) = incipient_kind
            end if

            return
         end if
      end do
   end subroutine check_critical_jump

   subroutine near_critical(nc, np, X, near_crit, l, i, j)
      integer, intent(in) :: nc !! Number of components
      integer, intent(in) :: np !! Number of main phases
      real(pr), intent(in) :: X(nc*np+np+2) !! Vector of variables
      logical, intent(out) :: near_crit
      !! If the system is near a critical point
      integer, intent(out) :: l !! Index of the near-critical phase
      integer, intent(out) :: i !! Index of the component with the maximum lnK
      integer, intent(out) :: j !! Index of the component with the minimum lnK


      integer :: lb, ub

      real(pr) :: cf, lnK(nc)

      near_crit = .false.

      do l=1,np
         lb = (l-1)*nc + 1
         ub = l*nc
         lnK = X(lb:ub)

         i = maxloc(lnK, dim=1)
         j = minloc(lnK, dim=1)
         near_crit= (maxval(exp(lnK)))/minval(exp(lnK)) - 1 < 0.2 .or. maxval(abs(lnK)) < 0.06_pr
         near_crit= (maxval(exp(lnK)))/minval(exp(lnK)) - 1 < 0.2 .or. maxval(abs(lnK)) < 0.1
         near_crit = maxval(lnK) - minval(lnK) < 0.5 .or. maxval(abs(lnK)) < 0.06
         if (near_crit) then
            return
         end if
      end do
   end subroutine near_critical
end module yaeos__equilibria_boundaries_auxiliar
