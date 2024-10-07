module yaeos__equilibria_saturation_points
   use yaeos__constants, only: pr
   use yaeos__models, only: ArModel
   use yaeos__equilibria_equilibrium_state, only: EquilibriumState
   use yaeos__equilibria_auxiliar, only: k_wilson
   use ieee_arithmetic, only: ieee_is_nan, ieee_is_finite

   implicit none

   real(pr) :: tol = 1e-9_pr
   integer :: max_iterations = 2000
   integer :: iters_first_step = 100
   real(pr) :: step_tol = 0.1_pr

   class(ArModel), pointer, private :: hidden_model
   real(pr), private, allocatable :: hidden_z(:)
   character(len=14), private :: hidden_kind

   real(pr), private :: Vz, Vy

contains

   type(EquilibriumState) function saturation_pressure(model, n, t, kind, p0, y0, max_iters)
      !! Saturation pressure calculation function.
      !!
      !! Calculates the saturation pressure of a multicomponent mixture with
      !! a given molar composition `n`.
      !! It is possible to calculate:
      !!
      !! - Bubble point: `kind="bubble"`
      !! - Dew point: `kind="dew"`
      !! - Liquid-Liquid point: `kind="liquid-liquid"`
      use stdlib_optval, only: optval
      class(ArModel), target, intent(in) :: model
      real(pr), intent(in) :: n(:) !! Composition vector [moles / molar fraction]
      real(pr), intent(in) :: t !! Temperature [K]
      character(len=*), intent(in) :: kind !! [bubble|dew|liquid-liquid]
      real(pr), optional, intent(in) :: p0 !! Initial pressure [bar]
      real(pr), optional, intent(in) :: y0(:) !! Initial composition
      integer, optional, intent(in) :: max_iters !! Maximum number of iterations

      real(pr) :: p, vy, vz

      real(pr) :: k(size(n)), y(size(n)), z(size(n)), lnk(size(n))
      real(pr) :: lnfug_y(size(n)), dlnphi_dp_y(size(n))
      real(pr) :: lnfug_z(size(n)), dlnphi_dp_z(size(n))

      character(len=50) :: incipient
      character(len=50) :: main

      real(pr) :: f, step
      integer :: its, iterations, i

      ! =======================================================================
      ! Handle arguments
      ! -----------------------------------------------------------------------
      z = n/sum(n)
      if (present (p0)) then
         p = p0
      else
         call model%pressure(z, T, 10._pr, P=P)
      end if

      if (present(y0)) then
         y = y0
      else
         y = z * k_wilson(model, T, P)
      end if
      iterations = optval(max_iters, max_iterations)

      select case(kind)
       case("bubble")
         k = y/z
         incipient = "vapor"
         main = "liquid"
       case("dew")
         k = z/y
         incipient = "liquid"
         main = "vapor"
       case("liquid-liquid")
         k = y/z
         incipient = "liquid"
         main = "liquid"
      end select

      where (z == 0)
         k = 0
      end where
      ! ========================================================================

      ! ========================================================================
      !  Solve point
      ! ------------------------------------------------------------------------
      do its=1, iters_first_step
         y = k*z
         call model%lnphi_pt(y, P, T, vy, incipient, lnPhi=lnfug_y, dlnphidp=dlnphi_dp_y)
         call model%lnphi_pt(z, P, T, vz, main, lnPhi=lnfug_z, dlnphidp=dlnphi_dp_z)

         k = exp(lnfug_z - lnfug_y)

         if (all(k < 1e-9_pr) .or. all(abs(k-1) < tol)) exit


         f = sum(z*k) - 1
         step = f/sum(z * k * (dlnphi_dp_z - dlnphi_dp_y))

         do while (P - step < 0 .or. abs(step) > 0.1*P)
            step = step/2
         end do

         p = p - step
         if (abs(step) < tol .and. abs(f) < tol) exit

      end do
      ! ========================================================================
      if (its > iters_first_step) then
      fsolve: block
         use yaeos__math_continuation, only: full_newton
         real(pr) :: X(size(n)+2)
         real(pr) :: S, dS, dXdS(size(n)+2)
         real(pr) :: F(size(n)+2), dF(size(n)+2, size(n)+2), dFdS(size(n)+2)
         integer :: ns

         ns = size(n)+1
         K = k_wilson(model, T, P)
         if (kind == "dew") K =1/K
         X = log([K, T, P])
         S = X(ns)

         hidden_kind = kind
         hidden_model => model
         hidden_z = z

         its = 0
         call full_newton(foo, its, X, ns, S, dS, dXdS, 1, max_iterations, F, dF, dFdS, tol=1.e-5_pr)
         K = exp(X(:size(n)))
         P = exp(X(size(n)+2))
         y = K*z
      end block fsolve
      end if
      
      select case(kind)
       case("bubble")
         saturation_pressure = EquilibriumState(kind="bubble", &
            iters=its, y=y, x=z, vx=vz, vy=vy, t=t, p=p, beta=0._pr&
            )
       case("dew")
         saturation_pressure = EquilibriumState(kind="dew", &
            iters=its, x=y, y=z, vy=vz, vx=vy, t=t, p=p, beta=1._pr&
            )
       case("liquid-liquid")
         saturation_pressure = EquilibriumState(kind="liquid-liquid", &
            iters=its, y=y, x=z, vx=vz, vy=vy, t=t, p=p, beta=0._pr&
            )
      end select
   end function saturation_pressure

   type(EquilibriumState) function saturation_temperature(model, n, p, kind, t0, y0, max_iters)
      !! Saturation temperature calculation function.
      !!
      !! Calculates the saturation pressure of a multicomponent mixture with
      !! a given molar composition `n`.
      !! It is possible to calculate:
      !!
      !! - Bubble point: `kind="bubble"`
      !! - Dew point: `kind="dew"`
      !! - Liquid-Liquid point: `kind="liquid-liquid"`
      use stdlib_optval, only: optval
      use yaeos__math_continuation, only: full_newton
      class(ArModel), target, intent(in) :: model
      real(pr), intent(in) :: n(:) !! Composition vector [moles / molar fraction]
      real(pr), intent(in) :: p !! Pressure [bar]
      character(len=*), intent(in) :: kind !! [bubble|dew|liquid-liquid]
      real(pr), optional, intent(in) :: t0 !! Initial temperature [K]
      real(pr), optional, intent(in) :: y0(:) !! Initial composition
      integer, optional, intent(in) :: max_iters !! Maximum number of iterations

      real(pr) :: t, vy, vz

      real(pr) :: k(size(n)), y(size(n)), z(size(n)), lnk(size(n))
      real(pr) :: lnfug_y(size(n)), dlnphi_dt_y(size(n))
      real(pr) :: lnfug_z(size(n)), dlnphi_dt_z(size(n))

      character(len=50) :: incipient
      character(len=50) :: main

      real(pr) :: f, step
      integer :: its, iterations

      logical :: is_incipient(size(n))

      ! =======================================================================
      ! Handle arguments
      ! -----------------------------------------------------------------------
      is_incipient = .true.
      z = n/sum(n)
      if (present (t0)) then
         t = t0
      else
         t = 150._pr
      end if

      if (present(y0)) then
         y = y0
      else
         y = z * k_wilson(model, T, P)
      end if
      iterations = optval(max_iters, max_iterations)

      select case(kind)
       case("bubble")
         k = y/z
         incipient = "vapor"
         main = "liquid"
       case("dew")
         k = z/y
         incipient = "liquid"
         main = "vapor"
       case("liquid-liquid")
         k = y/z
         incipient = "liquid"
         main = "liquid"
      end select

      where (z == 0)
         k = 0
      end where

      where (y == 0)
         is_incipient = .false.
      end where

      ! ========================================================================
      !  Solve point
      ! ------------------------------------------------------------------------
      do its=1, iters_first_step
         y = k*z
         where (.not. is_incipient)
            y = 0
         endwhere

         call model%lnphi_pt(y, P, T, vy, incipient, lnPhi=lnfug_y, dlnphidt=dlnphi_dt_y)
         call model%lnphi_pt(z, P, T, vz, main, lnPhi=lnfug_z, dlnphidt=dlnphi_dt_z)

         k = exp(lnfug_z - lnfug_y)
         f = sum(z*k) - 1
         step = f/sum(z * k * (dlnphi_dt_z - dlnphi_dt_y))

         if (.not. ieee_is_finite(step) .or. ieee_is_nan(step)) exit

         do while (abs(step) > 0.5*T .or. T - step < 0)
            if (isnan(step)) step = 10
            step = step/2
         end do

         t = t - step

         if (abs(step) < tol .and. abs(f) < tol) exit
      end do
      ! ========================================================================
      if (its > iters_first_step) then
      fsolve: block
         real(pr) :: X(size(n)+2)
         real(pr) :: S, dS, dXdS(size(n)+2)
         real(pr) :: F(size(n)+2), dF(size(n)+2, size(n)+2), dFdS(size(n)+2)
         integer :: ns

         ns = size(n)+2
         K = k_wilson(model, T, P)
         if (kind == "dew") K =1/K
         X = log([K, T, P])
         S = X(ns)

         hidden_kind = kind
         hidden_model => model
         hidden_z = z

         call full_newton(foo, its, X, ns, S, dS, dXdS, 1, max_iterations, F, dF, dFdS, tol=tol)
         K = exp(X(:size(n)))
         T = exp(X(size(n)+1))
         y = K*z
      end block fsolve
      end if

      select case(kind)
       case("bubble")
         saturation_temperature = EquilibriumState(kind="bubble", &
            iters=its, y=y, x=z, vx=vz, vy=vy, t=t, p=p, beta=0._pr&
            )
       case("dew")
         saturation_temperature = EquilibriumState(kind="dew", &
            iters=its, x=y, y=z, vy=vz, vx=vy, t=t, p=p, beta=1._pr&
            )
       case("liquid-liquid")
         saturation_temperature = EquilibriumState(kind="liquid-liquid", &
            iters=its, y=y, x=z, vx=vz, vy=vy, t=t, p=p, beta=0._pr&
            )
      end select


   end function saturation_temperature

   subroutine foo(X, ns, S, F, dF, dFdS)
      !! Function that needs to be solved at each envelope point
      real(pr), intent(in) :: X(:)
      integer, intent(in) :: ns
      real(pr), intent(in) :: S

      real(pr), intent(out) :: F(:)
      real(pr), intent(out) :: dF(:, :)
      real(pr), intent(out) :: dFdS(:)

      character(len=14) :: kind_z, kind_y

      real(pr) :: y(size(X)-2)
      real(pr) :: lnPhi_z(size(X)-2), lnPhi_y(size(X)-2)
      real(pr) :: dlnphi_dt_z(size(X)-2), dlnphi_dt_y(size(X)-2)
      real(pr) :: dlnphi_dp_z(size(X)-2), dlnphi_dp_y(size(X)-2)
      real(pr) :: dlnphi_dn_z(size(X)-2, size(X)-2), dlnphi_dn_y(size(X)-2, size(X)-2)

      real(pr) :: T, P, K(size(X)-2)

      real(pr) :: z(size(X)-2)

      integer :: i, j, nc

      nc = size(X)-2

      F = 0
      dF = 0

      K = exp(X(:nc))
      T = exp(X(nc+1))
      P = exp(X(nc+2))

      z = hidden_z
      y = K*z

      select case(hidden_kind)
       case ("bubble")
         kind_z = "liquid"
         kind_y = "vapor"
       case ("dew")
         kind_z = "vapor"
         kind_y = "liquid"
       case ("liquid-liquid")
         kind_z = "liquid"
         kind_y = "liquid"
       case default
         kind_z = "stable"
         kind_y = "stable"
      end select

      call hidden_model%lnphi_pt(&
         z, P, T, V=Vz, root_type=kind_z, &
         lnPhi=lnphi_z, dlnPhidt=dlnphi_dt_z, &
         dlnPhidp=dlnphi_dp_z, dlnphidn=dlnphi_dn_z &
         )
      call hidden_model%lnphi_pt(&
         y, P, T, V=Vy, root_type=kind_y, &
         lnPhi=lnphi_y, dlnPhidt=dlnphi_dt_y, &
         dlnPhidp=dlnphi_dp_y, dlnphidn=dlnphi_dn_y &
         )

      F(:nc) = X(:nc) + lnPhi_y - lnPhi_z
      F(nc + 1) = sum(y - z)
      F(nc + 2) = X(ns) - S

      ! Jacobian Matrix
      do j=1,nc
         df(:nc, j) = dlnphi_dn_y(:, j) * y(j)
         df(j, j) = dF(j, j) + 1
      end do

      df(:nc, nc + 1) = T * (dlnphi_dt_y - dlnphi_dt_z)
      df(:nc, nc + 2) = P * (dlnphi_dp_y - dlnphi_dp_z)

      df(nc + 1, :nc) = y

      df(nc + 2, :) = 0
      df(nc + 2, ns) = 1

      dFdS = 0
      dFdS(nc+2) = -1
   end subroutine foo
end module yaeos__equilibria_saturation_points
