module yaeos__equilibria_flash
   use yaeos__constants, only: pr
   use yaeos__models, only: ArModel
   use yaeos__equilibria_equilibria_state, only: EquilibriaState
   use yaeos__phase_equilibria_rachford_rice, only: betato01, betalimits, rachford_rice, solve_rr
   use yaeos__phase_equilibria_auxiliar, only: k_wilson
   implicit none

contains

   type(EquilibriaState) function flash(model, z, t, v_spec, p_spec, k0, iters)
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
      class(ArModel), intent(in) :: model !! Thermodynamic model
      real(pr), intent(in) :: z(:) !! Global composition (molar fractions)
      real(pr), intent(in) :: t !! Temperature [K]
      real(pr), optional, intent(in) :: v_spec !! Specified Volume [L/mol]
      real(pr), optional, intent(in) :: p_spec !! Specified Pressure [bar]
      real(pr), optional, intent(in) :: k0(:) !! Initial K factors (y/x)
      integer, optional, intent(out) :: iters !! Number of iterations

      logical :: stopflash

      ! Results from flash calculation
      real(pr), dimension(size(z)) :: x  ! composition of liquid (molar fractions)
      real(pr), dimension(size(z)) :: y  ! composition of vapour (molar fractions)
      real(pr) :: rho_x            ! density of liquid (moles/L)
      real(pr) :: rho_y            ! density of vapour (moles/L)
      real(pr) :: beta             ! total fraction of vapour (molar base)

      ! Intermediate variables during calculation process
      real(pr) :: p, v
      real(pr), dimension(size(z)) :: lnfug_y, lnfug_x
      real(pr), dimension(size(z)) :: K, dK, lnK, dKold, lnKold
      real(pr), dimension(size(z), size(z)) :: dlnphidn

      real(pr) :: g0, g1  ! function g valuated at beta=0 and 1, based on Wilson K factors
      real(pr) :: bmin, bmax, Vy, Vx

      real(pr) :: bx, savek(size(z)), log_k2(size(z))
      real(pr) :: DPVl, dpvv, dVydVl, h, pl, pold, pold2, pv, step, stepv
      real(pr) :: told, told2

      character(len=2) :: spec !! Flash specification [PT | VT]


      ! ========================================================================
      ! Starting steps
      ! ------------------------------------------------------------------------
      if (present(v_spec) .and. present(p_spec)) then
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
      do while (maxval(abs(dK)) > 1.d-6)
         iters = iters + 1

         call solve_rr(z, K, beta, bmin, bmax)

         y = z * K / (1 + beta*(K - 1._pr))
         x = y/K

         ! Calculate fugacities for each kind of specification
         select case (spec)
          case("TV")
            ! find Vy,Vx (vV and vL) from V balance and P equality equations
            call tv_loop_solve_pressures(model, T, V, beta, x, y, Vx, Vy, P)
            call model%lnphi_tp(y, T, P, V=Vy, root_type="stable", lnphip=lnfug_y)
            call model%lnphi_tp(x, T, P, V=Vx, root_type="liquid", lnphip=lnfug_x)
          case("TP")
            call model%lnphi_tp(y, T, P, V=Vy, root_type="stable", lnphip=lnfug_y)
            call model%lnphi_tp(x, T, P, V=Vx, root_type="liquid", lnphip=lnfug_x)
         end select

         dKold = dK
         lnKold = lnK

         lnK = lnfug_x - lnfug_y
         dK = lnK - lnKold

         K = exp(lnK)

         if (iters > 10 .and. abs(sum(dK + dKold)) < 0.05) then
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
            return
         end if
      end do

      ! ========================================================================
      ! Format results
      ! ------------------------------------------------------------------------
      if (spec == 'TP') v = beta*Vy + (1 - beta)*Vx

      if (maxval(K) < 1.001 .and. minval(K) > 0.999) then ! trivial solution
         flash%kind = "failed"
         P = -1.0
         flash%x = x/x
         flash%y = y/y
         flash%iters = iters
         flash%p = p
         flash%t = t
         return
      end if

      flash%kind = "split"
      flash%iters = iters
      flash%p = p
      flash%t = t

      flash%x = x
      flash%y = y
      flash%vx = Vx
      flash%vy = vy
      flash%beta = beta
   end function flash

   subroutine tv_loop_solve_pressures(model, T, V, beta, x, y, vx, vy, P)
      !! Solve pressure equality between two phases at a given temperature,
      !! total volume, vapor molar fractions and compositions.
      use iso_fortran_env, only: error_unit

      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: V !! Total volume [L/mol]
      real(pr), intent(in) :: beta !! Molar fraction of light-phase
      real(pr), intent(in) :: x(:) !! Molar fractions of heavy-phase
      real(pr), intent(in) :: y(:) !! Molar fractions of light-phase
      real(pr), intent(in out) :: Vx !! Heavy-phase molar volume [L/mol]
      real(pr), intent(in out) :: Vy !! Light-Phase molar volume [L/mol]
      real(pr), intent(out) :: P !! Pressure [bar]

      real(pr) :: Bx !! Liquid phase covolume
      real(pr) :: dVydVx !! Derivative of Vy wrt Vx

      ! Pressure equality newton functions
      real(pr) :: h !! Pressure equality
      real(pr) :: dh  !! dh/
      real(pr) :: stepv

      real(pr) :: dPxdV, dPydV
      real(pr) :: Px, Py

      integer :: its

      dVydVx = -(1 - beta)/beta
      Bx = model%get_v0(x, 0.1_pr, T)

      ! First evaluation will be with Vx = 1.5*Bx
      if (Vx < Bx) Vx = 1.625_pr*Bx

      call model%pressure(x, Vx, T, Px, dpdv=dPxdV)

      do while (Px < 0 .or. dPxdV >= 0)
         Vx = Vx - 0.2*(Vx - Bx)
         call model%pressure(x, Vx, T, Px, dpdv=dPxdV)
      end do

      Vy = (V - (1 - beta)*Vx)/beta

      h = 1.0
      its = 0
      do while (abs(h) > 1.d-4)
         ! Newton for solving P equality, with Vx as independent variable
         its = its + 1

         call model%pressure(x, Vx, T, Px, dpdv=dPxdV)
         call model%pressure(y, Vy, T, Py, dpdv=dPydV)

         h = Py - Px
         dh = -dPydV * dVydVx - dPxdV
         stepv = -h/dh

         if (its >= 10) stepv = stepv/2

         Vx = Vx + stepv

         do while (Vx < 1.001*Bx)
            stepv = stepv/2
            Vx = Vx - stepv
         end do

         Vy = (v - (1 - beta)*Vx)/beta

         if (its >= 100) then
            write (error_unit, *) "WARN(FLASH_VT): volume convergence problems", Px, Py
            P = -1.0
            return
         end if
      end do

      call model%pressure(x, Vx, T, Px)
      call model%pressure(y, Vy, T, Py)
      P = (Px + Py) * 0.5_pr
   end subroutine tv_loop_solve_pressures
end module yaeos__equilibria_flash
