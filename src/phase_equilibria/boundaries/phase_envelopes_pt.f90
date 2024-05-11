module yaeos__phase_equilibria_boundaries_phase_envelopes_pt
   !! Phase boundaries line on the \(PT\) plane calculation procedures.
   use yaeos_constants, only: pr
   use yaeos_models, only: ArModel
   use yaeos_equilibria_equilibria_state, only: EquilibriaState
   use yaeos_thermoprops, only: fugacity_tp
   use yaeos__math_continuation, only: continuation
   implicit none

   type :: CriticalPoint
      !! Critical point
      real(pr) :: T !! Temperature [K]
      real(pr) :: P !! Pressure [bar]
   end type CriticalPoint

   type :: PTEnvel2
      !! Two-phase isopleth.
      !! Phase boundary line of a fluid at constant composition.
      type(EquilibriaState), allocatable :: points(:)
      !! Each point through the line.
      type(CriticalPoint), allocatable :: cps(:)
      !! Critical points found along the line.
   end type PTEnvel2

contains

   function pt_envelope_2ph(&
      model, z, init_state, points, iterations &
      ) result(envelopes)
      !! PT two-phase envelope calculation procedure.
      !!
      !! Phase envelope calculation using the continuation method.
      use stdlib_optval, only: optval
      class(ArModel), intent(in) :: model !! Thermodyanmic model
      real(pr), intent(in) :: z(:) !! Vector of molar fractions
      type(EquilibriaState), intent(in) :: init_state !! Initial state
      integer, optional, intent(in) :: points !! Maxmimum number of points
      integer, optional, intent(in) :: iterations !! Point solver maxmimum iterations
      type(PTEnvel2) :: envelopes

      integer :: nc
      integer :: ns

      integer :: max_points
      integer :: its, max_iterations
      integer :: point

      real(pr), target :: X(size(z) + 2)
      real(pr), pointer :: lnK(:)
      real(pr), pointer :: lnT
      real(pr), pointer :: lnP
      real(pr), allocatable :: XS(:, :)
      real(pr) :: Xc(size(z) + 2)

      integer :: i, mc

      max_points = optval(points, 500)
      max_iterations = optval(iterations, 100)

      nc = size(z)
      lnK => X(:nc)
      lnT => X(nc + 1)
      lnP => X(nc + 2)

      lnK = log(init_state%y/init_state%x)
      lnT = log(init_state%T)
      lnP = log(init_state%P)

      XS = continuation(&
         foo, X, ns0=nc+1, S0=log(init_state%T), &
         dS0=0.1_pr, max_points=max_points, solver_tol=1.e-8_pr, &
         update_specification=update_specification &
         )

      allocate(envelopes%points(0), envelopes%cps(0))
      do i=2, size(XS, dim=1)
         if (all(XS(i, :) == 0._pr)) exit

         X = XS(i, :)

         envelopes%points = [&
            envelopes%points, &
            EquilibriaState(&
            x=z, Vx=0._pr, y=exp(lnK)*z, Vy=0, &
            T=exp(lnT), P=exp(lnP), beta=0._pr, iters=0) &
            ]

         if (all(XS(i - 1, :nc) * XS(i, :nc) < 0)) then
            mc = maxloc(abs(XS(i, :nc) - XS(i-1, :nc)), dim=1)
            Xc = &
               (XS(i, :) - XS(i-1, :))/(XS(i, mc) - XS(i-1,mc)) &
               * (XS(i, mc)) + (XS(i-1, :))
            envelopes%cps = [envelopes%cps, CriticalPoint(T=exp(lnT), P=exp(lnP))]
         end if
      end do

   contains
      subroutine foo(X, ns, S, F, dF, dFdS)
         !! Function that needs to be solved at each envelope point
         real(pr), intent(in) :: X(:)
         integer, intent(in) :: ns
         real(pr), intent(in) :: S

         real(pr), intent(out) :: F(:)
         real(pr), intent(out) :: dF(:, :)
         real(pr), intent(out) :: dFdS(:)

         real(pr) :: y(nc)
         real(pr) :: Vz, Vy, lnphip_z(nc), lnphip_y(nc)
         real(pr) :: dlnphi_dt_z(nc), dlnphi_dt_y(nc)
         real(pr) :: dlnphi_dp_z(nc), dlnphi_dp_y(nc)
         real(pr) :: dlnphi_dn_z(nc, nc), dlnphi_dn_y(nc, nc)

         real(pr) :: T, P, K(nc)

         integer :: i, j

         F = 0
         dF = 0

         K = exp(X(:nc))
         T = exp(X(nc+1))
         P = exp(X(nc+2))

         y = K*z

         call fugacity_tp(&
            model, z, T, P, V=Vz, root_type="stable", &
            lnphip=lnphip_z, dlnPhidt=dlnphi_dt_z, &
            dlnPhidp=dlnphi_dp_z, dlnphidn=dlnphi_dn_z &
            )
         call fugacity_tp(&
            model, y, T, P, V=Vy, root_type="stable", &
            lnphip=lnphip_y, dlnPhidt=dlnphi_dt_y, &
            dlnPhidp=dlnphi_dp_y, dlnphidn=dlnphi_dn_y &
            )

         F(:nc) = X(:nc) + lnphip_y - lnphip_z
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

      subroutine update_specification(X, ns, S, dS, dXdS, iterations)
         !! Update the specification during continuation.
         real(pr), intent(in out) :: X(:)
         !! Vector of variables \([\lnK_i \dots , lnT, lnP]\)
         integer, intent(in out) :: ns
         !! Number of specified variable in the vector
         real(pr), intent(in out) :: S
         !! Variable specification value
         real(pr), intent(in out) :: dS
         !! Step in specification
         real(pr), intent(in out) :: dXdS(:)
         !! Variation of variables with respect to specification
         integer, intent(in) :: iterations
         !! Iterations used in the solver

         real(pr) :: P, dP
         real(pr) :: T, dT
         real(pr) :: Xnew(nc+2)

         if (maxval(abs(X(:nc))) < 0.1_pr) then
            ns = maxloc(abs(dXdS(:nc)), dim=1)
         else
            ns = maxloc(abs(dXdS), dim=1)
         end if

         dS = dXdS(ns) * dS
         dXdS = dXdS/dXdS(ns)

         dS = sign(1.2_pr, dS) * minval([ &
            max(abs(sqrt(X(ns))/10), 0.1_pr), &
            abs(dS)*3/iterations &
            ] &
            )

         dS = sign(1._pr, dS) * maxval([abs(dS), 0.1_pr])
      end subroutine update_specification
   end function pt_envelope_2ph
end module yaeos__phase_equilibria_boundaries_phase_envelopes_pt
