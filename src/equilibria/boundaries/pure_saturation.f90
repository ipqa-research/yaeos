module yaeos__equilibria_boundaries_pure_saturation
   use yaeos__constants, only: pr
   use yaeos__models_ar, only: ArModel, size
   use yaeos__math_linalg, only: solve_system
   use yaeos__math_continuation, only: &
      continuation, continuation_solver, continuation_stopper
   use linear_interpolation_module, only: linear_interp_1d
   implicit none


   type :: PurePsat
      real(pr), allocatable :: T(:) !! Temperature [K]
      real(pr), allocatable :: P(:) !! Pressure [Pa]
      real(pr), allocatable :: Vx(:) !! Molar volume [L/mol] in the liquid phase
      real(pr), allocatable :: Vy(:) !! Molar volume [L/mol] in the vapor phase
      type(linear_interp_1d), private :: interpolator_get_T
      type(linear_interp_1d), private :: interpolator_get_P
   contains
      procedure :: get_T => get_T
      procedure :: get_P => get_P
   end type PurePsat

contains

   function pure_saturation_line(model, component, minP, minT) result(pt)
      !! # Pure saturation line
      !!
      !! Saturation pressures and temperatures for a pure component.
      !!
      !! ## Description
      !! This function calculates the saturation line for a pure component.
      !! Starting from the pure component critical point, the function traces
      !! the saturation line using the continuation method.
      !! The function returns a `PurePsat` object with the saturation
      !! temperatures and pressures. The object also contains interpolators
      !! to get the saturation temperature for a given pressure and vice versa.
      !!
      ! ========================================================================
      use yaeos__auxiliar, only: optval
      class(ArModel), intent(in) :: model !! Thermodyanmic model
      integer, intent(in) :: component !! Component index to calculate the line
      real(pr), intent(in) :: minP !! Minimum pressure [bar]
      real(pr), intent(in) :: minT !! Minimum temperature [K]
      type(PurePsat) :: pt
      ! ------------------------------------------------------------------------

      real(pr) :: X(4) !! Variables [lnVx, lnVy, lnP, lnT]
      real(pr) :: z(size(model))
      real(pr) :: Vc
      integer :: i
      integer :: ns

      real(pr) :: Tc, Pc
      real(pr) :: Vx, Vy, T, P

      real(pr) :: dXdS(4), dS, S, dFdS(4)
      real(pr) :: F(4), dF(4,4), dX(4)
      integer :: its, nc
      integer :: points

      integer :: iT, iP

      iP = 3
      iT = 4

      nc = size(model)
      Tc = model%components%Tc(component)
      Pc = model%components%Pc(component)

      z = 0
      z(component) = 1
      call model%volume(z, P=Pc, T=Tc, V=Vc, root_type="vapor")

      Vx = Vc*0.995
      Vy = Vc*1.005

      X = [log(Vx), log(Vy), log(Pc), log(Tc)]

      ns = 1
      S = log(0.95)
      dS = -0.001
      allocate(pt%T(0), pt%P(0), pt%Vx(0), pt%Vy(0))

      ! ========================================================================
      ! Trace the line using the continuation method.
      ! ------------------------------------------------------------------------
      T = Tc
      P = Pc
      points = 0
      do while(T > minT .and. P > minP .and. .not. isnan(T))
         call solve_point(model, component, nc, X, ns, S, F, dF, dFdS, its)
         dXdS = solve_system(dF, -dFdS)
         ns = maxloc(abs(dXdS(3:4)), dim=1) + 2
         dS = dXdS(ns)*dS
         dXdS = dXdS/dXdS(ns)

         do while (exp(X(4)) - exp(X(4) + dXdS(4)*dS) < 3 .and. ((Tc - T) > 10 .or. (Pc - P) > 2))
            dS = dS*1.5
         end do

         ds = sign(max(abs(dS),0.01_pr), dS)

         Vx = exp(X(1))
         Vy = exp(X(2))
         P  = exp(X(3))
         T  = exp(X(4))

         if (isnan(T)) then
            exit
         else
            pt%T = [pt%T, T]
            pt%P = [pt%P, P]
            pt%Vx = [pt%Vx, Vx]
            pt%Vy = [pt%Vy, Vy]
            points = points + 1
         end if

         dX = dXdS*dS
         do while(abs(exp(X(iT)) - exp(X(iT) + dX(iT))) > 10)
            dS = dS/2
            dX = dXdS*dS
         end do

         do while(abs(exp(X(iP)) - exp(X(iP) + dX(iP))) > 5)
            dS = dS/2
            dX = dXdS*dS
         end do

         do while(abs((exp(X(2)) - exp(X(2) + dX(2))) / exp(X(2))) > 0.1)
            dS = dS/2
            dX = dXdS*dS
         end do

         X = X + dX
         S = X(ns)
      end do

      ! Save interpolators to obtain particular values. The interpolator needs
      ! monothonic increasing values in x, so we need to reverse the arrays.
      pt%P  = pt%P(points:1:-1)
      pt%T  = pt%T(points:1:-1)
      pt%Vx = pt%Vx(points:1:-1)
      pt%Vy = pt%Vy(points:1:-1)

      call pt%interpolator_get_T%initialize(pt%P, pt%T, i)
      call pt%interpolator_get_P%initialize(pt%T, pt%P, i)
   end function pure_saturation_line

   subroutine solve_point(model, ncomp, nc, X, ns, S, F, dF, dFdS, its)
      !! # Solve point
      !!
      !! Solve a saturation point for a pure component.
      !!
      !! ## Description
      !! The set of equations to solve is:
      !!
      !! \[
      !! \begin{align*}
      !! f_1 &= \ln f_{z}(V_z, T) - \ln f_{y}(V_y, T) \\
      !! f_2 &= \ln \left( \frac{P_z}{P_y} \right) \\
      !! f_3 &= \ln P_z - \ln P \\
      !! f_4 &= g(X, ns)
      !! \end{align*}
      !! \]
      !!
      !! Where \(f_4\) is an specification function defined as:
      !!
      !! \[
      !! g(X, ns) = \left\{
      !! \begin{array}{lr}
      !! \ln \left( \frac{V_z}{V_y} \right) - S & \text{if } ns = 1 \text{ or } ns = 2 \\
      !! X(ns) - S & \text{otherwise}
      !! \end{array}
      !! \right\}
      !! \]
      !!
      !! The vector of variables \(X\) is equal to
      !! \([ \ln V_z, \ln V_y, \ln P, \ln T ]\).

      class(ArModel), intent(in) :: model
      !! Thermodynamic model
      integer, intent(in) :: ncomp
      !! Component index
      integer, intent(in) :: nc
      !! Total number of components
      real(pr), intent(in out) :: X(4)
      !! Variables \([ln V_z, lnV_y, lnP, lnT]\)
      integer, intent(in) :: ns
      !! Variable index to solve. If the
      real(pr), intent(in) :: S
      !! Variable value specified to solve
      real(pr), intent(out) :: F(4)
      !! Function
      real(pr), intent(out) :: dF(4, 4)
      !! Jacobian
      real(pr), intent(out) :: dFdS(4)
      !! Derivative of the function with respect to S
      integer, intent(out) :: its
      !! Number of iterations

      real(pr) :: z(nc)
      real(pr) :: lnfug_z(nc), lnfug_y(nc)
      real(pr) :: dlnfdv_z(nc), dlnfdv_y(nc)
      real(pr) :: dlnfdt_z(nc), dlnfdt_y(nc)
      real(pr) :: dPdTz, dPdTy
      real(pr) :: dPdVz, dPdVy
      real(pr) :: Vz, Vy

      real(pr) :: T
      real(pr) :: Pz, Py
      real(pr) :: dX(4), B
      real(pr) :: Xnew(4)

      integer :: i

      i = ncomp

      dX = 1
      F = 1
      z = 0
      z(i) = 1
      B = model%get_v0(z, 1._pr, 150._pr)

      its = 0
      do while((maxval(abs(dX)) > 1e-7 .and. maxval(abs(F)) > 1e-7))
         its = its+1
         call isofugacity(X, F, dF, dFdS)

         if (any(isnan(F))) exit
         dX = solve_system(dF, -F)
         Xnew = X + dX

         do while(abs(exp(Xnew(4)) - exp(X(4))) > 1)
            dX = dX/4
            Xnew = X + dX
         end do

         X = Xnew
      end do

   contains
      subroutine isofugacity(X, F, dF, dFdS)
         real(pr), intent(inout) :: X(4)
         real(pr), intent(out) :: F(4)
         real(pr), intent(out) :: dF(4,4)
         real(pr), intent(out) :: dFdS(4)

         F = 0
         dF = 0

         Vz = exp(X(1))
         Vy = exp(X(2))
         !lnP = X(3)
         T = exp(X(4))

         call model%lnfug_vt(z, V=Vz, T=T, P=Pz, lnf=lnfug_z, dlnfdV=dlnfdv_z, dlnfdT=dlnfdT_z, dPdV=dPdVz, dPdT=dPdTz)
         call model%lnfug_vt(z, V=Vy, T=T, P=Py, lnf=lnfug_y, dlnfdV=dlnfdv_y, dlnfdT=dlnfdT_y, dPdV=dPdVy, dPdT=dPdTy)

         F(1) = lnfug_z(i) - lnfug_y(i)
         F(2) = log(Pz/Py)
         F(3) = X(3) - log(Pz)

         if (ns == 1 .or. ns == 2) then
            F(4) = log(Vz/Vy) - S! X(ns) - S
         else
            F(4) = X(ns) - S
         end if

         dF = 0
         dF(1, 1) = Vz * dlnfdv_z(i)
         dF(1, 2) = -Vy * dlnfdv_y(i)
         dF(1, 3) = 0
         dF(1, 4) = T * (dlnfdT_z(i) - dlnfdT_y(i))

         dF(2, 1) = Vz/Pz * dPdVz
         dF(2, 2) = -Vy/Py * dPdVy
         dF(2, 4) = T * (dPdTz/Pz - dPdTy/Py)

         dF(3, 1) = -Vz/Pz * dPdVz
         dF(3, 2) = 0
         dF(3, 3) = 1
         dF(3, 4) = -T/Pz * dPdTz

         if (ns == 1 .or. ns == 2) then
            dF(4, 1) = 1
            dF(4, 2) = -1
         else
            dF(4, ns) = 1
         end if

         dFdS = 0
         dFdS(4) = -1
      end subroutine isofugacity
   end subroutine solve_point

   real(pr) function get_T(pt, P) result(T)
      !! # Get temperature
      !!
      !! Get the saturation temperature for a given pressure.
      !!
      !! ## Description
      !! This function returns the saturation temperature for a given pressure.
      !! The function uses an interpolator to get the required value.
      !!
      !! ## Examples
      !! ```fortran
      !! T = pt%get_T(P)
      !! ```
      class(PurePsat), intent(in out) :: pt
      real(pr), intent(in) :: P
      call pt%interpolator_get_T%evaluate(P, T)
   end function get_T

   real(pr) function get_P(pt, T) result(P)
      !! # Get pressure
      !!
      !! Get the saturation pressure for a given temperature.
      !!
      !! ## Description
      !! This function returns the saturation pressure for a given temperature.
      !! The function uses an interpolator to get the required value.
      !!
      !! ## Examples
      !! ```fortran
      !! P = pt%get_P(T)
      !! ```
      class(PurePsat), intent(in out) :: pt
      real(pr), intent(in) :: T
      call pt%interpolator_get_P%evaluate(T, P)
   end function get_P
end module yaeos__equilibria_boundaries_pure_saturation
