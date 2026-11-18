module yaeos__models_solvers
   !! # Model-specific numerical solvers
   !! Specialized numerical methods used by EOS implementations.
   !!
   !! # Description
   !! The principal routine is `volume_michelsen`. It solves in reduced
   !! covolume \(\zeta=B/V\), bounded between zero and one, rather than
   !! iterating directly over an unbounded volume. Models must implement
   !! `get_v0` to provide the covolume or lower-volume estimate.
   !!
   !! # Examples
   !! ```fortran
   !! call volume_michelsen(model, n, P, T, V, root_type="stable")
   !! ```
   !!
   !! # References
   !! Michelsen, M. L., and J. M. Mollerup, *Thermodynamic Models:
   !! Fundamentals & Computational Aspects*, 2nd ed., Tie-Line Publications, 2007.
   use yaeos__constants, only: pr, R, solving_volume
   use yaeos__models_ar, only: ArModel
   implicit none

contains

   subroutine volume_michelsen(eos, n, P, T, V, root_type, max_iters, V0)
      !! Volume solver at a given pressure.
      !!
      !! Obtain the volume using the method described by Michelsen and Mollerup.
      !! While \(P(V, T)\) can be obtained with a simple Newton method, a better
      !! approach is solving \(P(B/V, T)\) where \(B\) is the EoS covolume.
      !! This method is easier to solve because:
      !! \[
      !!    V(P, T) \in [0, \infty)
      !! \]
      !! and
      !! \[
      !!    \frac{B}{V}(P, T) \in [0, 1]
      !! \]
      !!
      !! At chapter 3 page 94 of Michelsen and Møllerup's book a more complete
      !! explanation can be seen
      !!
      !! # Notes
      !! Pass `root_type` explicitly (`"liquid"`, `"vapor"`, or `"stable"`).
      !! Although optional in the interface, the current implementation uses
      !! the argument directly in root selection. When `V0` is present only the
      !! initial \(\zeta = B/V_0\) is set and no root search is performed, so
      !! `V` is not assigned.
      use iso_fortran_env, only: error_unit
      use yaeos__auxiliar, only: optval
      class(ArModel), intent(in) :: eos !! EOS providing pressure and residual Helmholtz derivatives.
      real(pr), intent(in) ::  n(:) !! Mixture moles [mol]
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(out) :: V !! Volume [L]
      character(len=*), optional, intent(in) :: root_type !! Type of root ["vapor" | "liquid" | "stable"]
      integer, optional, intent(in) :: max_iters !! Maximum number of iterations (default 1000; currently unused, `solve_point` caps at 500)

      real(pr), optional, intent(in) :: V0 !! Initial volume guess [L]; see Notes

      character(len=10) :: root


      real(pr) :: totn !! Total number of moles.
      real(pr) :: B !! Covolume from `get_v0` [L]
      real(pr) :: ZETMIN
      real(pr) :: ZETA !! Reduced covolume \(\zeta = B/V\)
      real(pr) :: ZETMAX
      real(pr) :: pcalc
      real(pr) :: AT
      real(pr) :: AVAP
      real(pr) :: VVAP

      integer :: iter !! Iteration count.
      integer :: maximum_iterations

      maximum_iterations = optval(max_iters, 1000)
      root = optval(root_type, "stable")

      TOTN = sum(n)
      B = eos%get_v0(n, p, t)
      ITER = 0

      ! Limits
      ZETMIN = 0._pr
      ZETMAX = 1._pr
      ZETMAX = 1._pr - 0.01*T/(10000*B) ! improvement for cases with heavy components

      if (present(V0)) then
         zeta = B/V0
      else
         select case(root_type)
          case("liquid")
            ZETA = 0.5_pr
            call solve_point(eos, n, P, T, V, Pcalc, ZETA, ZETMIN, ZETMAX, AT, iter)
          case("vapor","stable")
            ZETA = min(0.5_pr, B*P/(TOTN*R*T))
            call solve_point(eos, n, P, T, V, Pcalc, ZETA, ZETMIN, ZETMAX, AT, iter)

            if (root_type == "stable") then
               ! Run first for vapor and then for liquid
               VVAP = V
               AVAP = AT
               ZETA = 0.5_pr
               ZETMAX = 1._pr
               ZETMAX = 1.D0 - 0.01*T/(10000*B) ! improvement for cases with heavy components
               call solve_point(eos, n, P, T, V, Pcalc, ZETA, ZETMIN, ZETMAX, AT, iter)
               if (AT .gt. AVAP) V = VVAP
            end if
          case default
            write(error_unit, *) "ERROR [VCALC]: Wrong specification", root_type
            error stop 1
         end select
      end if
   end subroutine volume_michelsen

   subroutine solve_point(eos, n, P, T, V, Pcalc, ZETA, ZETMIN, ZETMAX, AT, iter)
      !! # solve_point
      !! Correct one pressure root in reduced-covolume coordinates.
      !!
      !! # Description
      !! Updates \(\zeta=B/V\) with a bracketed Newton step, maintaining
      !! `ZETMIN` and `ZETMAX` as the current root bounds.
      class(ArModel), intent(in) :: eos !! EOS used for residual pressure evaluation.
      real(pr), intent(in) :: n(:) !! Component mole amounts.
      real(pr), intent(in) :: P !! Objective pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(out) :: V !! Obtained volume [L]
      real(pr), intent(out) :: Pcalc !! Calculated pressure at V [bar]
      real(pr), intent(in out) :: ZETA !! Current reduced covolume iterate, updated in place.
      real(pr), intent(inout) :: ZETMIN !! Lower bound for the reduced covolume root.
      real(pr), intent(inout) :: ZETMAX !! Upper bound for the reduced covolume root.
      real(pr), intent(out) :: AT !! \((A^r + PV)/(RT) - n\ln V\), a residual-Gibbs-like quantity used to select the stable root.
      integer, intent(out) :: iter !! Iterations used by the corrector.

      real(pr) :: del
      real(pr) :: der
      real(pr) :: B
      real(pr) :: totn !! Total number of moles.
      real(pr) :: Ar
      real(pr) :: ArV
      real(pr) :: ArV2 !! Second-order intermediate quantity used by this procedure.

      iter = 0
      DEL = 1
      pcalc = 2*p
      B = eos%get_v0(n, p, t)
      totn = sum(n)
      do while(abs(DEL) > 1.e-12_pr .and. iter < 500)
         V = B/ZETA
         iter = iter + 1

         solving_volume = .true.
         call eos%residual_helmholtz(n, V, T, Ar=Ar, ArV=ArV, ArV2=ArV2)
         solving_volume = .false.

         Pcalc = TOTN*R*T/V - ArV

         if (Pcalc .gt. P) then
            ZETMAX = ZETA
         else
            ZETMIN = ZETA
         end if

         ! AT is something close to Gr(P,T)
         AT = (Ar + V*P)/(T*R) - TOTN*log(V)

         ! this is dPdrho/B
         DER = (ArV2*V**2 + TOTN*R*T)/B

         DEL = -(Pcalc - P)/DER
         if (iter > 10) then
            ZETA = ZETA + max(min(DEL, 0.1_pr), -.1_pr)/2.
         else
            ZETA = ZETA + max(min(DEL, 0.1_pr), -.1_pr)
         end if

         if (ZETA .gt. ZETMAX .or. ZETA .lt. ZETMIN) then
            ZETA = 0.5_pr*(ZETMAX + ZETMIN)
         end if
      end do
   end subroutine solve_point
end module yaeos__models_solvers
