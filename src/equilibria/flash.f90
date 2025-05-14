module yaeos__equilibria_flash
   use yaeos__constants, only: pr
   use yaeos__models, only: BaseModel, ArModel, GeModel
   use yaeos__equilibria_equilibrium_state, only: EquilibriumState
   use yaeos__equilibria_rachford_rice, only: betato01, betalimits, rachford_rice, solve_rr
   use yaeos__equilibria_auxiliar, only: k_wilson
   use yaeos__solvers_pressure_equality, only: pressure_equality_V_beta_xy
   implicit none

contains

   type(EquilibriumState) function flash(model, z, t, v_spec, p_spec, k0, iters)
      !! Flash algorithm using sucessive substitutions.
      !!
      !! Available specifications:
      !!
      !! - TP (with T and P_spec variables)
      !! - TV (with T and V_spec variables)
      !!
      !! This algorithm assumes that the specified T and P/V correspond to
      !! vapor-liquid separation predicted by the provided model (0<beta<1) and
      !! solves the equilibria and mass-balance equations with a fixed-point
      !! method.
      use yaeos__auxiliar, only: optval
      class(BaseModel), intent(in) :: model !! Thermodynamic model
      real(pr), intent(in) :: z(:) !! Global composition (molar fractions)
      real(pr), intent(in) :: t !! Temperature [K]
      real(pr), optional, intent(in) :: v_spec !! Specified Volume [L/mol]
      real(pr), optional, intent(in) :: p_spec !! Specified Pressure [bar]
      real(pr), optional, intent(in) :: k0(:) !! Initial K factors (y/x)
      integer, optional, intent(out) :: iters !! Number of iterations

      ! Results from flash calculation
      real(pr), dimension(size(z)) :: x  ! composition of liquid (molar fractions)
      real(pr), dimension(size(z)) :: y  ! composition of vapour (molar fractions)
      real(pr) :: beta             ! total fraction of vapour (molar base)

      ! Intermediate variables during calculation process
      real(pr) :: P, V
      real(pr), dimension(size(z)) :: lnfug_y, lnfug_x
      real(pr), dimension(size(z)) :: K, dK, lnK, dKold, lnKold

      real(pr) :: g0, g1  ! function g valuated at beta=0 and 1, based on Wilson K factors
      real(pr) :: bmin, bmax, Vy, Vx

      character(len=2) :: spec !! Flash specification [PT | VT]


      ! ========================================================================
      ! Starting steps
      ! ------------------------------------------------------------------------
      if (present(V_spec) .and. present(P_spec)) then
         write (*, *) "ERROR: Can't specify pressure and volume in Flash"
         return
      else if (present(p_spec)) then
         spec = "TP"
         p = p_spec
      else if (present(v_spec)) then
         spec = "TV"
         v = v_spec
      end if

      if (spec == 'TV') then
         Vx = 0.0
         if (.not. present(k0)) then
            ! the EoS one-phase pressure will be used to estimate Wilson K factors
            call model%pressure(z, v_spec, t, p=p)
            if (P < 0) P = 1.0
         end if
      end if

      if (present(K0)) then
         K = K0
      else
         K = k_wilson(model, t, p)
      end if

      ! Get K values that assure that beta is between 0 and 1
      call betato01(z, K)

      ! now we must have  g0>0 and g1<0 and therefore 0<beta<1 (M&M page 252)
      call betalimits(z, K, bmin, bmax)
      beta = (bmin + bmax)/2  ! first guess for beta

      lnK = log(K)

      ! ========================================================================
      ! Solve with successive substitutions
      ! ------------------------------------------------------------------------
      dK = 1.0
      iters = 0
      do while (maxval(abs(dK)) > 1.e-6_pr)
         iters = iters + 1

         call betato01(z, K)
         call solve_rr(z, K, beta, bmin, bmax)

         y = z * K / (1 + beta*(K - 1._pr))
         x = y/K

         ! Calculate fugacities for each kind of specification

         select type (model)
          class is (GeModel)
            if (present(v_spec) .or. present(p_spec)) then
               error stop "Flash: GeModel can only spec T"
            end if
            call model%ln_activity_coefficient(y, T, lngamma=lnfug_y)
            call model%ln_activity_coefficient(x, T, lngamma=lnfug_x)
          class is (ArModel)
            select case (spec)
             case("TV")
               ! find Vy,Vx (vV and vL) from V balance and P equality equations
               call pressure_equality_V_beta_xy(model, T, V, beta, x, y, Vx, Vy, P)
               call model%lnphi_pt(y, P, T, V=Vy, root_type="stable", lnPhi=lnfug_y)
               call model%lnphi_pt(x, P, T, V=Vx, root_type="liquid", lnPhi=lnfug_x)
             case("TP")
               call model%lnphi_pt(y, P, T, V=Vy, root_type="stable", lnPhi=lnfug_y)
               call model%lnphi_pt(x, P, T, V=Vx, root_type="liquid", lnPhi=lnfug_x)
            end select
         end select

         dKold = dK
         lnKold = lnK

         lnK = lnfug_x - lnfug_y
         dK = lnK - lnKold

         K = exp(lnK)

         if (iters > 10 .and. abs(sum(dK + dKold)) < 0.05_pr) then
            ! oscilation behavior detected (27/06/15)
            lnK = (lnK + lnKold)/2
         end if

         ! Assure that beta is between the limits
         call betalimits(z, K, bmin, bmax)  ! 26/06/15
         if ((beta < bmin) .or. (bmax < beta)) then
            beta = (bmin + bmax)/2
         end if

         ! Step is too big, go back
         if (maxval(abs(dK)) > 1.10_pr) then
            ! 26/11/2014
            g0 = sum(z*K) - 1._pr
            g1 = 1._pr - sum(z/K)
            if (g0 < 0 .or. g1 > 0) then
               ! bring beta back to range, by touching K
               call betato01(z, K)
               call betalimits(z, K, bmin, bmax)
               beta = (bmin + bmax)/2  ! new guess for beta
            end if
         end if

         if (iters > 500) then
            p = -1
            exit
         end if
      end do

      ! ========================================================================
      ! Format results
      ! ------------------------------------------------------------------------
      if (spec == 'TP') V = beta*Vy + (1 - beta)*Vx

      if (maxval(K) < 1.001_pr .and. minval(K) > 0.999_pr .or. P < 0) then ! trivial solution
         flash%kind = "failed"
         P = -1.0
         flash%x = x/x
         flash%y = y/y
         flash%iters = iters
         flash%P = P
         flash%T = T
         return
      end if

      flash%kind = "split"
      flash%iters = iters
      flash%P = P
      flash%T = T

      flash%x = x
      flash%y = y
      flash%Vx = Vx
      flash%Vy = Vy
      flash%beta = beta
   end function flash
end module yaeos__equilibria_flash
