module yaeos__models_solvers
   !! # `models solvers`
   !! Set of different specialized solvers for different models
   !!
   !! # Description
   !! This module holds specialized solvers for different kind of applications
   !! and models.
   !!
   !! ## Volume solving
   !! This module holds the routine `volume_michelsen` which is a solver for
   !! volume that takes advantage over a simple newton on the function of
   !! pressure by solving the function of pressure over the covolume instead,
   !! which solution is limited in the range [0, 1]. This solver requires that
   !! the EoS uses the method `get_v0` to return the covolume.
   !!
   !! # Examples
   !!
   !! ```fortran
   !!  A basic code example
   !! ```
   !!
   !! # References
   !!
   use yaeos__constants, only: pr, R
   use yaeos__models_ar, only: ArModel
   implicit none

contains

   subroutine volume_michelsen(eos, n, P, T, V, root_type, max_iters)
      !! Volume solver at a given pressure.
      !!
      !! Obtain the volume using the method described by Michelsen and Møllerup.
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
      use iso_fortran_env, only: error_unit
      use stdlib_optval, only: optval
      class(ArModel), intent(in) :: eos
      real(pr), intent(in) ::  n(:) !! Mixture moles
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(out) :: V !! Volume [L]
      character(len=*), optional, intent(in) :: root_type !! Type of root ["vapor" | "liquid" | "stable"]
      integer, optional, intent(in) :: max_iters !! Maxiumum number of iterations, defaults to 100

      character(len=10) :: root

      real(pr) ::  Ar, ArV, ArV2

      real(pr) :: totn
      real(pr) :: B !! Covolume
      real(pr) :: ZETMIN, ZETA, ZETMAX
      real(pr) :: pcalc, AT, AVAP, VVAP

      integer :: iter, maximum_iterations

      maximum_iterations = optval(max_iters, 100)
      root = optval(root_type, "stable")

      TOTN = sum(n)
      B = eos%get_v0(n, p, t)
      ITER = 0

      ! Limits
      ZETMIN = 0._pr
      ZETMAX = 1._pr

      select case(root_type)
       case("liquid")
         ZETA = 0.5_pr
         call solve_point(P, V, Pcalc, AT, iter)
       case("vapor","stable")
         ZETA = min(0.5_pr, B*P/(TOTN*R*T))
         call solve_point(P, V, Pcalc, AT, iter)

         if (root_type == "stable") then
            ! Run first for vapor and then for liquid
            VVAP = V
            AVAP = AT
            ZETA = 0.5_pr
            ZETMAX = 1._pr
            call solve_point(P, V, Pcalc, AT, iter)
            if (AT .gt. AVAP) V = VVAP
         end if
       case default
         write(error_unit, *) "ERROR [VCALC]: Wrong specification"
         error stop 1
      end select
   contains
      subroutine solve_point(P, V, Pcalc, AT, iter)
         real(pr), intent(in) :: P !! Objective pressure [bar]
         real(pr), intent(out) :: V !! Obtained volume [L]
         real(pr), intent(out) :: Pcalc !! Calculated pressure at V [bar]
         real(pr), intent(out) :: AT !!
         integer, intent(out) :: iter

         real(pr) :: del, der

         iter = 0
         DEL = 1
         pcalc = 2*p
         do while(abs(DEL) > 1.e-10_pr .and. iter < maximum_iterations)
            V = B/ZETA
            iter = iter + 1
            call eos%residual_helmholtz(n, V, T, Ar=Ar, ArV=ArV, ArV2=ArV2)

            Pcalc = TOTN*R*T/V - ArV

            if (Pcalc .gt. P) then
               ZETMAX = ZETA
            else
               ZETMIN = ZETA
            end if

            ! AT is something close to Gr(P,T)
            AT = (Ar + V*P)/(T*R) - TOTN*log(V)

            DER = (ArV2*V**2 + TOTN*R*T)/B  ! this is dPdrho/B
            DEL = -(Pcalc - P)/DER
            ZETA = ZETA + max(min(DEL, 0.1_pr), -.1_pr)

            if (ZETA .gt. ZETMAX .or. ZETA .lt. ZETMIN) then
               ZETA = 0.5_pr*(ZETMAX + ZETMIN)
            end if
         end do

         if (iter >= maximum_iterations) write(error_unit, *) &
            "WARN: Volume solver exceeded maximum number of iterations"
      end subroutine solve_point
   end subroutine volume_michelsen
end module