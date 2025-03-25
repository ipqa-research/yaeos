module yaeos__equilibria_boundaries_phase_envelopes_pt3
   use yaeos__constants, only: pr, R
   use yaeos__equilibria_equilibrium_state, only: EquilibriumState
   use yaeos__models_ar, only: ArModel
   use yaeos__math, only: solve_system

   implicit none

   private

   public :: PTEnvel3
   public :: pt_envelope_3ph
   public :: solve_point
   public :: get_values_from_X

   type :: PTEnvel3
      integer, allocatable :: its(:) !! Number of needed iterations
      real(pr), allocatable :: beta(:) !! Mole fraction between phase x and phase y
      real(pr), allocatable :: x(:, :) !! Mole fraction of phase x
      real(pr), allocatable :: y(:, :) !! Mole fraction of phase x
      real(pr), allocatable :: w(:, :) !! Mole fraction of phase x
      real(pr), allocatable :: P(:) !! Pressures [bar]
      real(pr), allocatable :: T(:) !! Temperatures [K]
      integer, allocatable :: ns(:) !! Number of specified variable
      real(pr), allocatable :: S(:) !! Value of specification
   end type PTEnvel3

   real(pr), parameter :: lnK_min = 2.0_pr

contains

   type(PTEnvel3) function pt_envelope_3ph(&
      model, z, x0, y0, w0, beta0, P0, T0, ns0, dS0, &
      points &
      ) result(envelope)
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: z(:)
      real(pr), intent(in) :: x0(:)
      real(pr), intent(in) :: y0(:)
      real(pr), intent(in) :: w0(:)
      real(pr), intent(in) :: beta0
      real(pr), intent(in) :: P0
      real(pr), intent(in) :: T0
      integer, intent(in) :: ns0
      real(pr), intent(in) :: dS0
      integer, intent(in) :: points
      real(pr) :: kx(size(z))
      real(pr) :: ky(size(z))

      integer :: i
      integer :: nc

      real(pr) :: Xvars(size(z)*2 + 3), dX(size(z)*2 + 3)
      real(pr) :: F(size(z)*2+3), dF(size(z)*2 + 3, size(z)*2 + 3)

      integer :: ns !! Specified variable
      real(pr) :: S !! Specified value
      real(pr) :: dS !! Specified value step for next point extrapolation
      real(pr) :: dXdS(size(z)*2 + 3)

      real(pr) :: x(points, size(z)), xi(size(z))
      real(pr) :: y(points, size(z)), yi(size(z))
      real(pr) :: w(points, size(z)), wi(size(z))
      real(pr) :: beta(points)
      real(pr) :: P(points)
      real(pr) :: T(points)

      integer :: its, iterations(points)
      integer :: max_iterations

      max_iterations = 150

      nc = size(z)
      ns = ns0
      dS = dS0

      Kx = x0/w0
      Ky = y0/w0

      Xvars = [log(Kx), log(Ky), log(P0), log(T0), beta0]
      S = Xvars(ns)

      allocate(envelope%S(0), envelope%ns(0))

      ! Let the first point use more iterations
      call solve_point(model, z, ns, S, Xvars, F, dF, its, 1000)

      do i=1, points
         call solve_point(model, z, ns, S, Xvars, F, dF, its, max_iterations)
         if (any(isnan(F)) .or. any(isnan(Xvars)) .or. its >= max_iterations) exit

         envelope%ns = [envelope%ns, ns]
         envelope%S = [envelope%S, S]

         ! Set the values from the X vector.
         call get_values_from_X(z, Xvars, xi, yi, wi, P(i), T(i), beta(i))
         x(i, :) = xi
         y(i, :) = yi
         w(i, :) = wi

         iterations(i) = its
         call update_specification(its, Xvars, dF, dXdS, ns, dS)
         call detect_critical(Xvars, dXdS, ns, S, dS)

         dX = dXdS * dS
         Xvars = Xvars + dX
         S = Xvars(ns)
      end do

      i = i-1

      envelope%x = x(:i, :)
      envelope%y = y(:i, :)
      envelope%w = w(:i, :)
      envelope%P = P(:i)
      envelope%T = T(:i)
      envelope%beta = beta(:i)
      envelope%its = iterations(:i)
   end function pt_envelope_3ph

   subroutine get_values_from_X(z, Xvars, x, y, w, P, T, beta)
      real(pr), intent(in) :: z(:)
      real(pr), intent(in) :: Xvars(size(z)*2 + 3)
      real(pr), intent(out) :: x(size(z))
      real(pr), intent(out) :: y(size(z))
      real(pr), intent(out) :: w(size(z))
      real(pr), intent(out) :: P
      real(pr), intent(out) :: T
      real(pr), intent(out) :: beta

      integer :: nc
      real(pr) :: Kx((Size(Xvars)-3)/2), Ky((Size(Xvars)-3)/2)

      nc = (Size(Xvars)-3)/2

      Kx = exp(Xvars(1:nc))
      Ky = exp(Xvars(nc + 1:2*nc))
      P = exp(Xvars(2*nc + 1))
      T = exp(Xvars(2*nc + 2))
      beta = Xvars(2*nc + 3)

      w = z/(beta*Ky + (1 - beta)*Kx)
      x = w*Kx
      y = w*Ky
   end subroutine get_values_from_X

   subroutine update_specification(its, X, dF, dXdS, ns, dS)
      integer, intent(in) :: its
      real(pr), intent(in out) :: X(:)
      real(pr), intent(in out) :: dF(:, :)
      real(pr), intent(in out) :: dXdS(:)
      integer, intent(in out) :: ns
      real(pr), intent(in out) :: dS

      integer :: nc
      real(pr) :: dFdS(size(X))

      dFdS = 0
      dFdS(size(X)) = -1
      nc = (size(X) - 3)/2

      ns = maxloc(abs(dXdS), dim=1)

      dXdS = solve_system(dF, -dFdS)
      dS = dXdS(ns)*dS
      dXdS = dXdS/dXdS(ns)

      dS = dS * 3._pr/its

      do while(abs(dS/X(ns)) < 5e-2)
         dS = 2*dS
      end do
   end subroutine update_specification

   subroutine detect_critical(X, dXdS, ns, S, dS)
      !! # `detect_critical`
      !! Critical point detection
      !!
      !! # Description
      !! If the values of lnK (X[:nc]) change sign then a critical point
      !! Has passed, since for this to happen all variables should pass
      !! through zero. Near critical points (lnK < 0.05) points are harder
      !! to converge, so more steps in the extrapolation vector are made to
      !! jump over the critical point.
      !! If the critical point is detected then the kind of the point is
      !! changed and the point is saved using an interpolation knowing that
      !!
      !! \[
      !!   X_c = a * X + (1-a)*X_{new}
      !! \]
      !!
      !! With \(X_c\) is the variables at the critical point, \(X_{new}\)
      !! is the new initialization point of the method and \(a\) is the
      !! parameter to interpolate the values. This subroutine finds the
      !! value of  \(a\) to obtain \(X_c\).
      real(pr), intent(in out) :: X(:) !! Vector of variables
      real(pr), intent(in out) :: dXdS(:) !! Variation of variables wrt S
      integer, intent(in out) :: ns !! Number of specified variable
      real(pr), intent(in out) :: S !! Specification value
      real(pr), intent(in out) :: dS !! Step in specification
      real(pr) :: Xc(size(X)) !! Value at (near) critical point
      real(pr) :: a !! Parameter for interpolation

      real(pr) :: Xold(size(X)) !! Old value of X
      real(pr) :: Xnew(size(X)) !! Value of the next initialization

      logical :: found_critical

      integer :: nc

      integer :: first_set((size(X)-3)/2), second_set((size(X)-3)/2), idx((size(X)-3)/2)
      integer :: i, critical_set((size(X)-3)/2)
      real(pr) :: step(size(X))

      nc = (size(X)-3)/2
      first_set  = [(i, i=1, nc)]
      second_set = [(i, i=nc+1, 2*nc)]

      Xold = X
      Xnew = X + dXdS*dS

      found_critical = .false.

      do i=1,2
         select case(i)
          case(1)
            idx = first_set
          case(2)
            idx = second_set
         end select

         do while(maxval(abs(Xnew(idx))) < 0.3)
            Xnew = Xnew + dXdS*dS
         end do

         if (all(Xnew(idx) * Xold(idx) < 0)) then
            ! If two steps imply the crossing of a critical point, then
            ! make those two steps to avoid falling into it
            found_critical = .true.
            critical_set = idx
            Xnew = Xnew + dXdS*dS
         end if
      end do

      if (found_critical) then
         a = critical_interpol(Xnew, Xold, critical_set)
         Xc = a * Xold + (1-a)*Xnew

         ! Xnew = Xc
         ! do while(maxval(abs(Xnew(idx))) < 0.5)
         !    Xnew = Xnew + dXdS*dS
         ! end do
         X = X + (Xnew - Xold)
      end if
   end subroutine detect_critical

   real(pr)  function critical_interpol(Xnew, Xold, idx) result(a)
      !! # `critical_interpol`
      !! Critical point interpolation
      !!
      !! # Description
      !! This function calculates the parameter \(a\) to interpolate the
      !! values of the variables at the critical point. The interpolation
      !! is done using the equation:
      !!
      !! \[
      !!   0 = a*X_{old}(ns) + (1-a)*X_{new}(ns)
      !! \]
      !!
      !! Where \(X_{old}\) is the old value of the variables and \(X_{new}\)
      !! is the new value of the variables. The critical point is the point
      !! where the variables change sign, so the interpolation is done to
      !! find the value of the variables at the critical point.
      real(pr), intent(in) :: Xnew(:) !! New value of the variables
      real(pr), intent(in) :: Xold(:) !! Old value of the variables
      integer, intent(in) :: idx(:) !! Index of the variables to interpolate

      integer :: ncomp

      ! 0 = a*X(ns) + (1-a)*Xnew(ns) < Interpolation equation to get X(ns) = 0

      ! Find the component that is changing sign with the highest slope
      ncomp = maxloc(abs(Xold(idx) - Xnew(idx)), dim=1)
      a = -Xnew(ncomp)/(Xold(ncomp) - Xnew(ncomp))
   end function critical_interpol

   subroutine pt_F_three_phases(model, z, Xvars, ns, S, F, dF)
      !! Function to solve at each point of a three phase envelope.
      !!
      !! The vector of variables X corresponds to:
      !! \( X = [lnKx_i, lnKy_i lnP, lnT, \beta] \)
      !!
      !! While the equations are:
      !!
      !! \( F = [
      !!        lnKx_i - ln \phi_i(x, P, T) + ln \phi_i(w, P, T),
      !!        lnKy_i - ln \phi_i(y, P, T) + ln \phi_i(w, P, T),
      !!        \sum_{i=1}^N (w_i) - 1,
      !!        \sum_{i=1}^N (x_i - y_i),
      !!        X_{ns} - S
      !! ] \)
      use iso_fortran_env, only: error_unit
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: z(:)
      real(pr), intent(in)  :: Xvars(:) !! Vector of variables
      integer, intent(in)  :: ns !! Number of specification
      real(pr), intent(in)  :: S !! Specification value
      real(pr), intent(out) :: F(size(Xvars)) !! Vector of functions valuated
      real(pr), intent(out) :: df(size(Xvars), size(Xvars)) !! Jacobian matrix

      ! Xvars variables
      real(pr) :: Kx((Size(Xvars)-3)/2)
      real(pr) :: Ky((Size(Xvars)-3)/2)
      real(pr) :: P
      real(pr) :: T
      real(pr) :: beta

      ! Main phase 1 variables
      real(pr) :: Vx
      real(pr), dimension((Size(Xvars)-3)/2) :: x, lnphi_x, dlnphi_dt_x, dlnphi_dp_x
      real(pr), dimension((Size(Xvars)-3)/2, (Size(Xvars)-3)/2) :: dlnphi_dn_x

      ! Main phase 2 variables
      real(pr) :: Vy
      real(pr), dimension((Size(Xvars)-3)/2) :: y, lnphi_y, dlnphi_dt_y, dlnphi_dp_y
      real(pr), dimension((Size(Xvars)-3)/2, (Size(Xvars)-3)/2) :: dlnphi_dn_y

      ! Incipient phase variables
      real(pr) :: Vw
      real(pr), dimension((Size(Xvars)-3)/2) :: w, lnphi_w, dlnphi_dt_w, dlnphi_dp_w
      real(pr), dimension((Size(Xvars)-3)/2, (Size(Xvars)-3)/2) :: dlnphi_dn_w

      ! Derivative of w wrt beta
      real(pr) :: dwdb((Size(Xvars)-3)/2)

      real(pr) :: dwdKx((Size(Xvars)-3)/2), dxdKx((Size(Xvars)-3)/2), dydKx((Size(Xvars)-3)/2)
      real(pr) :: dwdKy((Size(Xvars)-3)/2), dxdKy((Size(Xvars)-3)/2), dydKy((Size(Xvars)-3)/2)

      integer :: i, j, nc

      nc = (size(Xvars)-3)/2

      Kx = exp(Xvars(1:nc))
      Ky = exp(Xvars(nc + 1:2*nc))
      P = exp(Xvars(2*nc + 1))
      T = exp(Xvars(2*nc + 2))
      beta = Xvars(2*nc + 3)

      w = z/(beta*Ky + (1 - beta)*Kx)
      x = w*Kx
      y = w*Ky

      call model%lnphi_pt(&
         x, P, T, V=Vx, root_type="stable", lnphi=lnphi_x, &
         dlnphidp=dlnphi_dp_x, dlnphidt=dlnphi_dt_x, dlnphidn=dlnphi_dn_x &
         )
      call model%lnphi_pt(&
         y, P, T, V=Vy, root_type="stable", lnphi=lnphi_y, &
         dlnphidp=dlnphi_dp_y, dlnphidt=dlnphi_dt_y, dlnphidn=dlnphi_dn_y &
         )
      call model%lnphi_pt(&
         w, P, T, V=Vw, root_type="stable", lnphi=lnphi_w, &
         dlnphidp=dlnphi_dp_w, dlnphidt=dlnphi_dt_w, dlnphidn=dlnphi_dn_w &
         )


      F = 0
      df = 0

      F(1:nc) = Xvars(1:nc) + lnphi_x - lnphi_w
      F(nc + 1:2*nc) = Xvars(nc + 1:2*nc) + lnphi_y - lnphi_w

      F(2*nc + 1) = sum(w) - 1
      F(2*nc + 2) = sum(x - y)
      F(2*nc + 3) = Xvars(ns) - S

      dwdb = z*(Kx - Ky)/((1 - beta)*Kx + beta*Ky)**2

      dwdKx = -z*(1 - beta)/(Ky*beta + (1 - beta)*Kx)**2
      dxdKx = Kx*dwdKx + w
      dydKx = Ky*dwdKx

      dwdKy = -z*(beta)/(Ky*beta + (1 - beta)*Kx)**2
      dxdKy = Kx*dwdKy
      dydKy = Ky*dwdKy + w

      do i = 1, nc
         do j = 1, nc
            df(i, j) = Kx(j)*(dlnphi_dn_x(i, j)*dxdKx(j) &
               - dlnphi_dn_w(i, j)*dwdKx(j))
            df(i + nc, j) = Kx(j)*(dlnphi_dn_y(i, j)*dydKx(j) &
               - dlnphi_dn_w(i, j)*dwdKx(j))

            df(i, j + nc) = Ky(j)*(dlnphi_dn_x(i, j)*dxdKy(j) &
               - dlnphi_dn_w(i, j)*dwdKy(j))
            df(i + nc, j + nc) = Ky(j)*(dlnphi_dn_y(i, j)*dydKy(j) &
               - dlnphi_dn_w(i, j)*dwdKy(j))
         end do

         ! dlnK_i/dlnK_i
         df(i, i) = df(i, i) + 1
         df(i + nc, i + nc) = df(i + nc, i + nc) + 1

         df(i, 2*nc + 3) = sum(Kx*dlnphi_dn_x(i, :)*dwdb - dlnphi_dn_w(i, :)*dwdb)
         df(i + nc, 2*nc + 3) = sum(Ky*dlnphi_dn_y(i, :)*dwdb - dlnphi_dn_w(i, :)*dwdb)

         df(2*nc + 1, i) = Kx(i)*dwdKx(i)
         df(2*nc + 1, i + nc) = Ky(i)*dwdKy(i)

         df(2*nc + 2, i) = Kx(i)*dxdKx(i) - Kx(i)*dydKx(i)
         df(2*nc + 2, i + nc) = Ky(i)*dxdKy(i) - Ky(i)*dydKy(i)
      end do

      ! Derivatives wrt P
      df(:nc, 2*nc + 1) = P*(dlnphi_dp_x - dlnphi_dp_w)
      df(nc + 1:2*nc, 2*nc + 1) = P*(dlnphi_dp_y - dlnphi_dp_w)

      ! Derivatives wrt T
      df(:nc, 2*nc + 2) = T*(dlnphi_dt_x - dlnphi_dt_w)
      df(nc + 1:2*nc, 2*nc + 2) = T*(dlnphi_dt_y - dlnphi_dt_w)

      ! Derivatives wrt beta
      df(2*nc + 1, 2*nc + 3) = sum(dwdb)
      df(2*nc + 2, 2*nc + 3) = sum(Kx*dwdb - Ky*dwdb)

      ! Derivatives wrt Xs
      df(2*nc + 3, :) = 0
      df(2*nc + 3, ns) = 1
   end subroutine pt_F_three_phases

   subroutine solve_point(model, z, ns, S, X, F, dF, its, maxits)
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: z(:)
      integer, intent(in) :: ns
      real(pr), intent(in) :: S
      real(pr), intent(in out) :: X(:)
      real(pr), intent(out) :: F(:)
      real(pr), intent(out) :: dF(:,:)
      integer, intent(in out) :: its
      integer, intent(in) :: maxits

      real(pr) :: dX(size(X))
      integer :: nc, i

      its = 0
      F = 1
      dX = 1
      nc = (size(X) - 3)/2


      do while((maxval(abs(F)) > 1e-7 .and. its < maxits))

         its = its + 1

         call pt_F_three_phases(model, z, X, ns, S, F, dF)

         dX = solve_system(dF, -F)

         do while((abs(dX(2*nc+1)/X(2*nc+1))) > 0.1)
            dX = dX/2
         end do

         do while((abs(dX(2*nc+2)/X(2*nc+2))) > 0.1)
            dX = dX/2
         end do

         do while(abs(dX(2*nc+3)) > 0.1)
            dX = dX/2
         end do

         X = X + dX

      end do
   end subroutine solve_point
end module yaeos__equilibria_boundaries_phase_envelopes_pt3
