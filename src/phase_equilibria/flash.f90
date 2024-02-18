module yaeos_equilibria_flash
   use yaeos_constants, only: pr
   use yaeos_models, only: ArModel
   use yaeos_equilibria_equilibria_state, only: EquilibriaState
   use yaeos_thermoprops, only: fugacity_vt, fugacity_tp, pressure
   implicit none

contains

   type(EquilibriaState) function flash(self, z, t, v_spec, p_spec, k0, iters)
      !! This algorithm assumes that the specified T and P correspond to
      !! vapor-liquid separation predicted by the provided model (0<beta<1)
      class(ArModel), intent(in) :: self !! Thermodynamic model
      real(pr), intent(in) :: z(:) !! Global composition (molar fractions)
      real(pr), intent(in) :: t !! Temperature [K]
      real(pr), optional, intent(in) :: v_spec !! Specified Volume [L/mol]
      real(pr), optional, intent(in) :: p_spec !! Specified Pressure [bar]
      real(pr), intent(in) :: k0(:) !! Initial K factors (y/x)
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

      real(pr) :: aux, bx, savek(size(z)), log_k2(size(z))
      real(pr) :: DPVl, dpvv, dVydVl, h, pl, pold, pold2, pv, step, stepv
      real(pr) :: told, told2

      character(len=2) :: spec !! Flash specification [PT | VT]


      ! ========================================================================
      ! Starting steps
      ! ------------------------------------------------------------------------
      if (present(v_spec) .and. present(p_spec)) then
         write (*, *) "ERROR: Can't specify pressure and volume in Flash"
      else if ( present(v_spec) ) then
         spec = "TV"
         v = v_spec
      else if ( present(p_spec) ) then
         spec = "TP"
         p = p_spec
      end if

      if (spec == 'TV' .or. spec == 'isoV') then
         Vx = 0.0
         ! if (FIRST) then  
            ! the EoS one-phase pressure will be used to estimate Wilson K factors
            call pressure(self, z, v_spec, t, p=p)
            if (P < 0) P = 1.0
         ! end if
      end if

      K = K0
      
      call betato01(z, K)  ! adapted 26/11/2014
      lnK = log(K)
      
      ! now we must have  g0>0 and g1<0 and therefore 0<beta<1 (M&M page 252)
      call betalimits(z, K, bmin, bmax)
      beta = (bmin + bmax)/2  ! first guess for beta
      ! ========================================================================

      ! Succesive sustitution loop starts here
      dK = 1.0
      iters = 0
      do while (maxval(abs(dK)) > 1.d-6)
         iters = iters + 1

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

         call solve_rr(z, K, beta, bmin, bmax)

         y = z * K / (1 + beta*(K - 1._pr))
         x = y/K

         ! new for TV Flash
         select case (spec)
         case("TV", "isoV")
            ! find Vy,Vx (vV and vL) from V balance and P equality equations
            ! TODO: Add TV specification
         case("TP")
            ! for TP Flash
            call fugacity_tp(self, y, T, P, V=Vy, root_type="stable", lnfug=lnfug_y)
            call fugacity_tp(self, x, T, P, V=Vx, root_type="liquid", lnfug=lnfug_x)
         end select

         dKold = dK
         lnKold = lnK

         lnK = lnfug_x - lnfug_y
         dK = lnK - lnKold
         
         aux = sum(dK + dKold)

         if (iters > 10 .and. abs(aux) < 0.05) then 
            ! oscilation behavior detected (27/06/15)
            lnK = (lnK + lnKold)/2
         end if

         K = exp(lnK)
         
         call betalimits(z, K, bmin, bmax)  ! 26/06/15
         if ((beta < bmin) .or. (bmax < beta)) then
            beta = (bmin + bmax)/2
         end if

         if (iters > 500) then
            p = -1
            return
         end if

      end do

      if (spec == 'TP') v = beta*Vy + (1 - beta)*Vx
      if (spec == 'TV' .or. spec == 'isoV') write (4, *) T, P, Pv
      if (spec == 'TV' .or. spec == 'isoV') P = Pv
      
      if (maxval(K) < 1.001 .and. minval(K) > 0.999) then ! trivial solution
         P = -1.0
         flash%x = x/x
         flash%y = y/y
         flash%iters = iters
         flash%p = p
         flash%t = t
         return
      end if

      flash%iters = iters
      flash%p = p
      flash%t = t

      if (Vy < Vx) then
         ! `y` phase is the heaviest phase
         flash%y = x
         flash%x = y
         flash%vy = Vx
         flash%vx = vy
      else
         flash%x = x
         flash%y = y
         flash%vx = Vx
         flash%vy = vy
      end if
   end function flash

   subroutine betato01(z, K)
      implicit none
      real(pr), intent(in) :: z(:) ! composition of the system
      real(pr) :: K(:) ! K factors (modified in this routine)
      real(pr) :: g0, g1  ! function g valuated at beta=0 and 1, based on K factors

      g1 = 1.0
      do while (g0 < 0 .or. g1 > 0)
         g0 = sum(z*K) - 1.D0
         g1 = 1.D0 - sum(z/K)
         if (g0 < 0) then
            K = 1.1*K  ! increased volatiliy will bring the solution from subcooled liquid into VLE
         else if (g1 > 0) then
            K = 0.9*K  ! decreased volatiliy will bring the solution from superheated vapor into VLE
         end if
      end do
   end subroutine betato01

   subroutine betalimits(z, K, bmin, bmax)
      real(pr), intent(in) :: z(:), K(:)  ! composition of the system and K factors
      real(pr), intent(out) :: bmin, bmax
      real(pr), dimension(size(z)) :: vmin, vmax
      integer :: i, in, ix

      in = 0
      ix = 0
      vmin = 0.d0
      ! max=1.001d0    ! modified  3/3/15 (not to generate false separations with beta 0.9999...)
      vmax = 1.00001d0 ! modified 28/6/15 (to prevent overshooting in the Newton for solving RR eq.)
      do i = 1, size(z)
         if (K(i)*z(i) > 1) then
            in = in + 1
            vmin(in) = (K(i)*z(i) - 1.d0)/(K(i) - 1.d0)
         else if (K(i) < z(i)) then
            ix = ix + 1
            vmax(ix) = (1.d0 - z(i))/(1.d0 - K(i))
         end if
      end do
      bmin = maxval(vmin)
      bmax = minval(vmax)
   end subroutine betalimits

   subroutine rachford_rice(z, K, beta, rr, drrdb)
      real(pr), intent(in) :: z(:)
      real(pr), intent(in) :: K(:)
      real(pr), intent(in) :: beta

      real(pr), intent(out) :: rr 
      real(pr), intent(out) :: drrdb

      real(pr) :: denom(size(z))
      
      denom = 1 + beta*(K - 1._pr)
      rr = sum(z*(K - 1._pr)/denom)
      drrdb = -sum(z*(K - 1._pr)**2/denom**2)
   end subroutine

   subroutine solve_rr(z, K, beta, beta_min, beta_max)
      !! Solve the Rachford-Rice Equation.
      real(pr), intent(in) :: z(:) !! Mole fractions vector
      real(pr), intent(in) :: K(:) !! K-factors
      real(pr), intent(out) :: beta_min !! 
      real(pr), intent(out) :: beta_max
      real(pr), intent(out) :: beta

      real(pr) :: g, dgdb
      real(pr) :: step

      g = 1.0
      step = 1.0

      call betalimits(z, k, beta_min, beta_max)

      do while (abs(g) > 1.d-5 .and. abs(step) > 1.d-10)
         call rachford_rice(z, k, beta, g, dgdb)
         step = -g/dgdb
         beta = beta + step
         do while ((beta < beta_min .or. beta_max < beta) .and. step > 1e-10) 
            step = step/2
            beta = beta - step
         end do
      end do
   end subroutine

   ! subroutine tv_loop
   !    dVydVl = -(1 - beta)/beta
   !    ! call Bcalc(n, x, T, Bx)
   !    ! TODO: Add this intiial volume

   !    if (Vx < Bx) Vx = 1.625*Bx  ! First evaluation will be with Vx = 1.5*Bx
   !    ! Pl = -1.0
   !    call zTVTERMO(n, 0, T, x, Vx, Pl, DPVl, lnfug_y, dlnphidp, dlnphitp, dlnphidn)  ! 26/06/15
   !    do while (Pl < 0 .or. DPVl >= 0)
   !       Vx = Vx - 0.2*(Vx - Bx)
   !       call zTVTERMO(n, 0, T, x, Vx, Pl, DPVl, lnfug_y, dlnphidp, dlnphitp, dlnphidn)
   !    end do
   !    Vy = (v - (1 - beta)*Vx)/beta
   !    h = 1.0
   !    iterv = 0

   !    stopflash = .false.
   !    do while (abs(h) > 1.d-4)  ! Newton for solving P equality, with Vx as independent variable
   !       iterv = iterv + 1
   !       if (iterv >= 100) then
   !          write (2, *) 'volume convergence problems'
   !          P = -1.0
   !          stopflash = .true.
   !          exit
   !       end if
   !       call zTVTERMO(n, 0, T, x, Vx, Pl, DPVl, lnfug_y, dlnphidp, dlnphitp, dlnphidn)
   !       call zTVTERMO(n, 0, T, y, Vy, Pv, DPVv, lnfug_y, dlnphidp, dlnphitp, dlnphidn)
   !       h = Pv - Pl
   !       dh = -DPVv*dVydVl - DPVl
   !       stepv = -h/dh
   !       if (iterv >= 10) stepv = stepv/2
   !       Vx = Vx + stepv
   !       do while (Vx < 1.001*Bx)
   !          stepv = stepv/2
   !          Vx = Vx - stepv
   !       end do
   !       Vy = (v - (1 - beta)*Vx)/beta
   !    end do
   !    if (stopflash .eqv. .true.) exit

   !    call zTVTERMO(n, 1, T, x, Vx, Pl, DPVl, PHILOGx, dlnphidp, dlnphitp, dlnphidn)
   !    call zTVTERMO(n, 1, T, y, Vy, Pv, DPVv, lnfug_y, dlnphidp, dlnphitp, dlnphidn)
   ! end subroutine
end module
