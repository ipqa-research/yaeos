module phase_equilibria
   use constants, only: pr
   use legacy_ar_models, only: zTVTERMO, termo, n => nc, omg => w, tc, pc
   implicit none

contains

   subroutine flash_pt(z, p, t, x, y, beta, its)
      real(pr), intent(in) :: z(:) !! Feed phase molar fractions
      real(pr), intent(in) :: p !! Pressure [bar]
      real(pr), intent(in) :: t !! Temperature [K]
      real(pr), intent(out) :: x(size(z)) !! Phase X molar fractions
      real(pr), intent(out) :: y(size(z)) !! Phase Y molar fractions
      real(pr), intent(out) :: beta !! Molar Y fraction
      real(pr), optional, intent(out) :: its !! Number of iterations

      real(pr) :: rho_x, rho_y
   end subroutine

   subroutine flash(spec, FIRST, z, t, p, v, x, y, rho_x, rho_y, beta, iter)
      ! Flash specification, eos id and  number of compounds in the system
      character(len=*), intent(in) :: spec !! Flash specification [PT | VT]
      logical, intent(in out) ::  FIRST
      logical :: stopflash

      ! composition of the system
      real*8, intent(in) :: z(:)

      ! Temperature and Pressure for the flash
      real*8, intent(in) :: t            ! Temperature for the flash (K)
      real*8 :: p            ! (bar) Pressure for the flash (TP) or resulting from (TV)
      real*8 :: v          ! (L/mol) Molar vol for the flash (TV) or resulting from (TP)

      ! Results from flash calculation
      real*8, dimension(size(z)), intent(out) :: x  ! composition of liquid (molar fractions)
      real*8, dimension(size(z)), intent(out) :: y  ! composition of vapour (molar fractions)
      real*8, intent(out) :: rho_x            ! density of liquid (moles/L)
      real*8, intent(out) :: rho_y            ! density of vapour (moles/L)
      real*8, intent(out) :: beta             ! total fraction of vapour (molar base)
      integer, intent(out) :: iter            ! number of iterations required to converge

      ! Intermediate variables during calculation process
      real*8, dimension(n) :: PHILOGy, PHILOGx, DLPHIT, DLPHIP
      real*8, dimension(n) :: KFACT, LOG_K, AUXK, var_K, denom, varKold, logKold
      real*8, dimension(n, n) :: FUGN
      real*8 :: g0, g1  ! function g valuated at beta=0 and 1, based on Wilson K factors
      real*8 :: g, dg, bmin, bmax, Vy, Vx

      ! real*8, dimension(nco, nco) :: Kij_or_K0, Tstar
      ! real*8, dimension(nco) :: saveK, LOG_K2
      real(8) :: aux, bx, savek(n), log_k2(n)
      integer :: MTYP

      real(8) :: dh, dpv, DPVl, dpvv, dVydVl, h, pl, pold, pold2, pv, step, stepv
      real(8) :: told, told2, bij(n, n)

      integer :: i, j, iterv, nco

      ! do i = 1, n
      !    do j = i, n
      !       bij(i, j) = (1 - lij(i, j))*(b(i) + b(j))/2
      !       bij(j, i) = bij(i, j)
      !    end do
      ! end do
      !
      !-----------------------------------------------------------
      ! This algorithm assumes that the specified T and P correspond to
      ! vapor-liquid separation predicted by the provided model (0<beta<1)

      if (spec == 'TV' .or. spec == 'isoV') then
         Vx = 0.0
         if (FIRST) then  ! the EoS one-phase pressure will be used to estimate Wilson K factors
            call zTVTERMO(n, 0, T, z, V, P, DPV, PHILOGy, DLPHIP, DLPHIT, FUGN)
            if (P < 0) P = 1.0
         end if
      end if
      AUXK = log(saveK(1:n))
      if (FIRST) then          !  use Wilson to initiate the first flash
         KFACT = (PC/P)*exp(5.373*(1 + omg)*(1 - TC/T))
         Pold2 = 0.d0
         Pold = 0.d0
         Told2 = 0.d0
         Told = 0.d0
      else !  for running the indirect "Tv flash" for comparisonn purposes
   !        else if(Pold2==0.d0.or.spec=='TV')then ! use the converged K's from the previous flash
         KFACT = saveK(1:n)
   !        else ! use extrapolation based on the last two points (not resolved yet for series of TV flashes)
   !            if(spec=='isoV')    LOG_K = AUXK + (AUXK - LOG_K2)
   !            if(spec=='TP'.and.P/=Pold)      LOG_K = AUXK + (AUXK - LOG_K2)*(P-Pold)/(Pold-Pold2)
   !            if(spec=='TP'.and.T/=Told)      LOG_K = AUXK + (AUXK - LOG_K2)*(T-Told)/(Told-Told2)
   !            KFACT = EXP(LOG_K)
      end if
      LOG_K2(:n) = AUXK
      Pold2 = Pold
      Pold = P
      Told2 = Told
      Told = T
      ! WRITE (2,3) (KFACT(i),i=1,N)
      call betato01(n, z, KFACT)  ! adapted 26/11/2014
      LOG_K = log(KFACT)
      ! now we must have  g0>0 and g1<0 and therefore 0<beta<1 (M&M page 252)
      call betalimits(n, z, KFACT, bmin, bmax)
      beta = (bmin + bmax)/2  ! first guess for beta

      ! Succesive sustitution loop starts here
      var_K = 1.0
      iter = 0
      do while (maxval(abs(var_K)) > 1.d-6)
         if (maxval(abs(var_K)) > 1.10) then  ! 26/11/2014
            g0 = sum(z*KFACT) - 1.D0
            g1 = 1.D0 - sum(z/KFACT)

            if (g0 < 0 .or. g1 > 0) then  ! bring beta back to range, by touching KFACT
               call betato01(n, z, KFACT)
               call betalimits(n, z, KFACT, bmin, bmax)
               beta = (bmin + bmax)/2  ! new guess for beta
            end if

         end if

         iter = iter + 1
         ! Newton starts here (Rachford-Rice)
         g = 1.0
         step = 1.0

         do while (abs(g) > 1.d-5 .and. abs(step) > 1.d-10)
            denom = 1 + beta*(KFACT - 1.D0)
            g = sum(z*(KFACT - 1.D0)/denom)
            dg = -sum(z*(KFACT - 1.D0)**2/denom**2)
            step = -g/dg
            beta = beta + step

            do while ((beta < bmin .or. bmax < beta) .and. step > 1e-10) ! much better (GUARANTED!) 3/3/15
               step = step/2
               beta = beta - step
            end do

         end do

         denom = 1 + beta*(KFACT - 1.D0)
         y = z*KFACT/denom
         x = y/KFACT

         ! new for TV Flash
         if (spec == 'TV' .or. spec == 'isoV') then     ! find Vy,Vx (vV and vL) from V balance and P equality equations
            dVydVl = -(1 - beta)/beta
            ! call Bcalc(n, x, T, Bx)
            ! TODO: Add this intiial volume

            if (Vx < Bx) Vx = 1.625*Bx  ! First evaluation will be with Vx = 1.5*Bx
            ! Pl = -1.0
            call zTVTERMO(n, 0, T, x, Vx, Pl, DPVl, PHILOGy, DLPHIP, DLPHIT, FUGN)  ! 26/06/15
            do while (Pl < 0 .or. DPVl >= 0)
               Vx = Vx - 0.2*(Vx - Bx)
               call zTVTERMO(n, 0, T, x, Vx, Pl, DPVl, PHILOGy, DLPHIP, DLPHIT, FUGN)
            end do
            Vy = (v - (1 - beta)*Vx)/beta
            h = 1.0
            iterv = 0

            stopflash = .false.
            do while (abs(h) > 1.d-4)  ! Newton for solving P equality, with Vx as independent variable
               iterv = iterv + 1
               if (iterv >= 100) then
                  write (2, *) 'volume convergence problems'
                  P = -1.0
                  stopflash = .true.
                  exit
               end if
               call zTVTERMO(n, 0, T, x, Vx, Pl, DPVl, PHILOGy, DLPHIP, DLPHIT, FUGN)
               call zTVTERMO(n, 0, T, y, Vy, Pv, DPVv, PHILOGy, DLPHIP, DLPHIT, FUGN)
               h = Pv - Pl
               dh = -DPVv*dVydVl - DPVl
               stepv = -h/dh
               if (iterv >= 10) stepv = stepv/2
               Vx = Vx + stepv
               do while (Vx < 1.001*Bx)
                  stepv = stepv/2
                  Vx = Vx - stepv
               end do
               Vy = (v - (1 - beta)*Vx)/beta
            end do
            if (stopflash .eqv. .true.) exit

            call zTVTERMO(n, 1, T, x, Vx, Pl, DPVl, PHILOGx, DLPHIP, DLPHIT, FUGN)
            call zTVTERMO(n, 1, T, y, Vy, Pv, DPVv, PHILOGy, DLPHIP, DLPHIT, FUGN)
         else  
            ! for TP Flash
            ! nc,MTYP,INDIC,T,P,rn,V,PHILOG,DLPHI
            MTYP = 0    ! -1   (with 0, generalized also fo LL and not only VL)
            call TERMO(n, MTYP, 1, T, P, y, Vy, PHILOGy, DLPHIP, DLPHIT, FUGN)
            MTYP = 1
            call TERMO(n, MTYP, 1, T, P, x, Vx, PHILOGx, DLPHIP, DLPHIT, FUGN)
         end if

         varKold = var_K
         logKold = LOG_K ! From previous iteration step
         var_K = PHILOGx - PHILOGy - LOG_K  ! variation in LOG_K = new - old
         LOG_K = PHILOGx - PHILOGy
         aux = sum(var_K + varKold)

         if (iter > 10 .and. abs(aux) < 0.05) then ! oscilation behavior detected (27/06/15)
            LOG_K = (LOG_K + logKold)/2
         end if

         KFACT = exp(LOG_K)
         call betalimits(n, z, KFACT, bmin, bmax)  ! 26/06/15

         if ((beta < bmin) .or. (bmax < beta)) then
            beta = (bmin + bmax)/2
         end if

         if (iter > 500) then
            p = -1
            return
         end if

      end do

      !  WRITE (2,4) (KFACT(i),i=1,N)
      rho_x = 1/Vx
      rho_y = 1/Vy
      if (spec == 'TP') v = beta*Vy + (1 - beta)*Vx
      if (spec == 'TV' .or. spec == 'isoV') write (4, *) T, P, Pv
      if (spec == 'TV' .or. spec == 'isoV') P = Pv
      FIRST = .false.

      if (maxval(KFACT) < 1.001 .and. minval(KFACT) > 0.999) then ! trivial solution
         P = -1.0
         return
         !go to 31
      end if
      saveK(1:n) = KFACT
   ! 3  format('KWilson ', 15E12.4)
   ! 4  format('KFinal  ', 15E12.4)
      !-----------------------------------------------------------

      ! print *, x  ! Estos print son los que "lee" tanto Fluids como Sur
      ! print *, y
      ! print *, rho_x
      ! print *, rho_y
      ! print *, beta
   end subroutine flash

   subroutine betato01(n, z, KFACT)
      implicit none
      integer, intent(in) :: n  ! number of compounds in the system
      real*8, dimension(n), intent(in) :: z ! composition of the system
      real*8, dimension(n) :: KFACT  ! K factors (modified in this routine)
      real*8 :: g0, g1  ! function g valuated at beta=0 and 1, based on K factors

      g1 = 1.0
      do while (g0 < 0 .or. g1 > 0)
         g0 = sum(z*KFACT) - 1.D0
         g1 = 1.D0 - sum(z/KFACT)
         if (g0 < 0) then
            KFACT = 1.1*KFACT  ! increased volatiliy will bring the solution from subcooled liquid into VLE
         else if (g1 > 0) then
            KFACT = 0.9*KFACT  ! decreased volatiliy will bring the solution from superheated vapor into VLE
         end if
      end do
   end subroutine betato01

   subroutine betalimits(n, z, KFACT, bmin, bmax)
      implicit none
      integer, intent(in) :: n  ! number of compounds in the system
      real*8, dimension(n), intent(in) :: z, KFACT  ! composition of the system and K factors
      real*8, intent(out) :: bmin, bmax
      real*8, dimension(n) :: vmin, vmax
      integer :: i, in, ix

      in = 0
      ix = 0
      vmin = 0.d0
      ! max=1.001d0    ! modified  3/3/15 (not to generate false separations with beta 0.9999...)
      vmax = 1.00001d0 ! modified 28/6/15 (to prevent overshooting in the Newton for solving RR eq.)
      do i = 1, n
         if (KFACT(i)*z(i) > 1) then
            in = in + 1
            vmin(in) = (KFACT(i)*z(i) - 1.d0)/(KFACT(i) - 1.d0)
         else if (KFACT(i) < z(i)) then
            ix = ix + 1
            vmax(ix) = (1.d0 - z(i))/(1.d0 - KFACT(i))
         end if
      end do
      bmin = maxval(vmin)
      bmax = minval(vmax)
   end subroutine betalimits
end module