module yaeos__equilibria_boundaries_phase_envelopes_px
   !! Phase boundaries line on the \(P\alpha\) plane calculation procedures.
   use yaeos__constants, only: pr
   use yaeos__models, only: ArModel
   use yaeos__equilibria_equilibrium_state, only: EquilibriumState
   use yaeos__math_continuation, only: &
      continuation, continuation_solver, continuation_stopper
   implicit none

   type :: CriticalPoint
      !! Critical point
      real(pr) :: alpha !! \(\alpha\)
      real(pr) :: P !! Pressure [bar]
   end type CriticalPoint

   type :: PXEnvel2
      !! Two-phase PX envelope.
      !! Phase boundary line of a fluid at constant temperature
      !! with variation in composition.
      real(pr), allocatable :: alpha(:) !! Second fluid molar fraction
      real(pr), allocatable :: z0(:) !! Original fluid composition
      real(pr), allocatable :: z_inj(:) !! Second fluid composition
      type(EquilibriumState), allocatable :: points(:)
      !! Each point through the line.
      type(CriticalPoint), allocatable :: cps(:)
      !! Critical points found along the line.
   end type PXEnvel2


   ! Private volumes of each phase to share between functions
   real(pr), private :: Vz !! Main phase volume [L/mol]
   real(pr), private :: Vy !! Incipient phase volume [L/mol]

contains

   function px_envelope_2ph(&
      model, z0, alpha0, z_injection, first_point, &
      points, iterations, delta_0, specified_variable_0, &
      solver, stop_conditions &
      ) result(envelopes)
      !! PX two-phase envelope calculation procedure.
      !!
      !! Phase envelope calculation using the continuation method.
      !! Defaults to solving the saturation temperature and continues with
      !! an increment in it. The variable to specify can be changed by modifying
      !! `specified_variable_0` with the corresponding variable number.
      ! ========================================================================
      use stdlib_optval, only: optval
      class(ArModel), intent(in) :: model
      !! Thermodyanmic model
      real(pr), intent(in) :: z0(:)
      !! Vector of molar fractions of the global composition (main phase)
      real(pr), intent(in) :: alpha0
      !! First point of \(alpha\)
      real(pr), intent(in) :: z_injection(:)
      !! Vector of molar fractions of the injection fluid
      type(EquilibriumState) :: first_point
      integer, optional, intent(in) :: points
      !! Maxmimum number of points, defaults to 500
      integer, optional, intent(in) :: iterations
      !! Point solver maximum iterations, defaults to 100
      real(pr), optional, intent(in) :: delta_0
      !! Initial extrapolation \(\Delta\)
      integer, optional, intent(in) :: specified_variable_0
      !! Position of specified variable, since the vector of variables is
      !! \(X = [lnK_i, \dots, lnP, \alpha]\) the values for specification
      !! will be \([1 \dots nc]\) for the equilibria constants, \(nc+1\) for
      !! \(lnP\) and \(nc + 2\) for \(\alpha\).
      procedure(continuation_solver), optional :: solver
      !! Specify solver for each point, defaults to a full newton procedure
      procedure(continuation_stopper), optional :: stop_conditions
      !! Function that returns true if the continuation method should stop
      type(PXEnvel2) :: envelopes
      ! ------------------------------------------------------------------------

      integer :: nc !! Number of components
      integer :: ns !! Number of specified variable
      real(pr) :: dS0 !! Initial specification step
      real(pr) :: S0 !! Initial specification value
      real(pr) :: z(size(z0)) !! Composition at some point

      integer :: max_points !! Maximum number of points
      integer :: max_iterations !! Maximum number of iterations

      real(pr) :: X(size(z) + 2), T
      real(pr), allocatable :: XS(:, :)

      character(len=14) :: kind

      ! ========================================================================
      ! Handle input
      ! ------------------------------------------------------------------------
      call get_z(alpha0, z0, z_injection, z)
      kind = first_point%kind
      nc = size(z)
      max_points = optval(points, 500)
      max_iterations = optval(iterations, 100)
      ns = optval(specified_variable_0, nc+2)
      dS0 = optval(delta_0, 0.1_pr)

      ! Correctly define the K-values based on the provided incipient point.
      select case(first_point%kind)
       case("bubble", "liquid-liquid")
         X(:nc) = log(first_point%y/z)
       case("dew")
         X(:nc) = log(first_point%x/z)
      end select

      T = first_point%T

      X(nc+1) = log(first_point%P)
      X(nc+2) = alpha0
      S0 = X(ns)
      ! print *, X, S0
      ! call exit
      allocate(envelopes%points(0), envelopes%cps(0), envelopes%alpha(0))

      test_numdiff: block
         real(pr) :: F(size(X)), df(size(X), size(X)), numdiff(size(X), size(X))
         real(pr) :: FdX(size(X)), dx(size(X)), dFdS(size(X))
         real(pr) :: FdX2(size(X))
         integer :: i
         integer :: loc(2)
         real(pr) :: maxerr
         exit test_numdiff

         do i=1,size(X)
            dx = 0
            dx(i) = 1.e-5_pr * X(i)
            call foo(X - dx, ns, S0, FdX, df, dFdS)
            call foo(X + dx, ns, S0, FdX2, df, dFdS)
            call foo(X, ns, S0, F, df, dFdS)
            numdiff(:, i) = (FdX2 - FdX)/(2*dx(i))
         end do

         loc = maxloc(abs(numdiff - df))
         maxerr = abs(&
            (numdiff(loc(1), loc(2)) - df(loc(1), loc(2))&
            )/numdiff(loc(1), loc(2)))
         if (maxerr > 0.01_pr) then
            print *, "ERROR: PXEnvel2 Numerical differentiation failed"
            loc = maxloc(abs(numdiff - df))
            print *, loc
            print *, df(loc(1), loc(2)), numdiff(loc(1), loc(2))
            error stop 1
         end if
      end block test_numdiff


      ! ========================================================================
      ! Trace the line using the continuation method.
      ! ------------------------------------------------------------------------
      XS = continuation(&
         foo, X, ns0=ns, S0=S0, &
         dS0=dS0, max_points=max_points, solver_tol=1.e-5_pr, &
         update_specification=update_spec, &
         solver=solver, stop=stop_conditions &
         )

   contains

      recursive subroutine foo(X, ns, S, F, dF, dFdS)
         !! Function that needs to be solved at each envelope point
         real(pr), intent(in) :: X(:)
         integer, intent(in) :: ns
         real(pr), intent(in) :: S

         real(pr), intent(out) :: F(:)
         real(pr), intent(out) :: dF(:, :)
         real(pr), intent(out) :: dFdS(:)

         character(len=14) :: kind_z, kind_y

         real(pr) :: y(nc)
         real(pr) :: lnphip_z(nc), lnphip_y(nc)
         real(pr) :: dlnphi_dt_z(nc), dlnphi_dt_y(nc)
         real(pr) :: dlnphi_dp_z(nc), dlnphi_dp_y(nc)
         real(pr) :: dlnphi_dn_z(nc, nc), dlnphi_dn_y(nc, nc)

         real(pr) :: P, K(nc), alpha, dzda(nc)

         integer :: i, j

         F = 0
         dF = 0

         K = exp(X(:nc))
         P = exp(X(nc+1))
         alpha = X(nc+2)

         call get_z(alpha, z0, z_injection, z, dzda)

         y = K*z

         select case(kind)
          case ("bubble")
            kind_z = "liquid"
            kind_y = "vapor"
          case ("dew")
            kind_z = "vapor"
            kind_y = "liquid"
          case default
            kind_z = "stable"
            kind_y = "stable"
         end select

         call model%lnphi_pt(&
            z, P=P, T=T, V=Vz, root_type=kind_z, &
            lnphi=lnphip_z, dlnPhidt=dlnphi_dt_z, &
            dlnPhidp=dlnphi_dp_z, dlnphidn=dlnphi_dn_z &
            )
         call model%lnphi_pt(&
            y, P=P, T=T, V=Vy, root_type=kind_y, &
            lnphi=lnphip_y, dlnPhidt=dlnphi_dt_y, &
            dlnPhidp=dlnphi_dp_y, dlnphidn=dlnphi_dn_y &
            )

         F(:nc) = X(:nc) + lnphip_y - lnphip_z
         F(nc + 1) = sum(y - z)
         F(nc + 2) = X(ns) - S

         ! Jacobian Matrix
         do i = 1, nc
            do j = 1, nc
               df(i, j) = y(j) * dlnphi_dn_y(i, j)
            end do
            df(i, i) = df(i, i) + 1
            
            df(i, nc + 2) = sum(K*dlnphi_dn_y(i, :)*dzda - dlnphi_dn_z(i, :)*dzda)
         end do

         df(:nc, nc + 1) = P*(dlnphi_dp_y - dlnphi_dp_z)
         df(nc + 1, :nc) = y
         df(nc + 1, nc + 2) = sum(dzda*(K - 1))

         df(nc + 2, :) = 0
         df(nc + 2, ns) = 1

         dFdS = 0
         dFdS(nc+2) = -1
      end subroutine foo

      subroutine update_spec(X, ns, S, dS, dXdS, step_iters)
         !! Update the specification during continuation.
         real(pr), intent(in out) :: X(:)
         !! Vector of variables \([lnK_i \dots , lnT, lnP]\)
         integer, intent(in out) :: ns
         !! Number of specified variable in the vector
         real(pr), intent(in out) :: S
         !! Variable specification value
         real(pr), intent(in out) :: dS
         !! Step in specification
         real(pr), intent(in out) :: dXdS(:)
         !! Variation of variables with respect to specification
         integer, intent(in) :: step_iters
         !! Iterations used in the solver

         real(pr) :: maxdS

         ! =====================================================================
         ! Update specification
         ! - Dont select T or P near critical points
         ! - Update dS wrt specification units
         ! - Set step
         ! ---------------------------------------------------------------------
         if (maxval(abs(X(:nc))) < 0.1_pr .and. abs(Vz - Vy) < 0.05) then
            ns = maxloc(abs(dXdS(:nc)), dim=1)
            maxdS = 0.01_pr
         else
            ns = maxloc(abs(dXdS), dim=1)
            maxdS = 0.5_pr
         end if

         dS = dXdS(ns) * dS
         dXdS = dXdS/dXdS(ns)
         

         dS = sign(1.0_pr, dS) * minval([ &
            max(sqrt(abs(X(ns))/10._pr), 0.1_pr), &
            abs(dS)*3/step_iters &
            ] &
            )

         ! do while (maxval(abs(dXdS(:nc)*dS)) < 0.01)
         !    dS = dS*1.1_pr
         ! end do

         call save_point(X, step_iters)
         call detect_critical(X, dXdS, ns, S, dS)
      end subroutine update_spec

      subroutine save_point(X, iters)
         !! Save the converged point
         real(pr), intent(in) :: X(:)
         integer, intent(in) :: iters
         type(EquilibriumState) :: point

         real(pr) :: y(nc), P, alpha

         P = exp(X(nc+1))
         alpha = X(nc+2)
         y = exp(X(:nc))*z

         select case(kind)
          case("bubble")
            point = EquilibriumState(&
               kind=kind, x=z, Vx=Vz, y=y, Vy=Vy, &
               T=T, P=P, beta=0._pr, iters=iters &
               )
          case("dew")
            point = EquilibriumState(&
               kind=kind, x=y, Vx=Vy, y=z, Vy=Vz, &
               T=T, P=P, beta=0._pr, iters=iters &
               )
          case default
            point = EquilibriumState(&
               kind="saturation", x=z, Vx=Vz, y=y, Vy=Vy, &
               T=T, P=P, beta=0._pr, iters=iters &
               )
         end select

         envelopes%alpha = [envelopes%alpha, alpha]
         envelopes%points = [envelopes%points, point]
      end subroutine save_point

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
         real(pr) :: Xc(nc+2) !! Value at (near) critical point
         real(pr) :: a !! Parameter for interpolation

         real(pr) :: Xold(size(X)) !! Old value of X
         real(pr) :: Xnew(size(X)) !! Value of the next initialization

         Xold = X

         do while (maxval(abs(X(:nc))) < 0.1_pr .and. abs(Vz - Vy) < 0.05_pr)
            ! If near a critical point, jump over it
            if (nc == 2) exit
            S = S + dS
            X = X + dXdS*dS
         end do

         Xnew = X + dXdS*dS

         if (all(Xold(:nc) * (Xnew(:nc)) < 0)) then
            

            select case(kind)
             case("dew")
               kind = "bubble"
             case("bubble")
               kind = "dew"
             case default
               kind = "liquid-liquid"
            end select

            ! 0 = a*X(ns) + (1-a)*Xnew(ns) Interpolation equation to get X(ns) = 0
            a = -Xnew(ns)/(X(ns) - Xnew(ns))
            Xc = a * X + (1-a)*Xnew

            envelopes%cps = [&
               envelopes%cps, &
               CriticalPoint(P=exp(Xc(nc+1)), alpha=Xc(nc+2)) &
               ]
            X = Xc + dXdS*dS
            
            if (nc == 2) then
               X = Xc
               dS = 0
            end if

         end if
      end subroutine detect_critical
   end function px_envelope_2ph

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

end module yaeos__equilibria_boundaries_phase_envelopes_px
