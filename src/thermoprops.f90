module yaeos_thermoprops
   !! Residual thermodyamic properties using residual Helmholtz model
   use yaeos_constants, only: R, pr
   use yaeos_models_ar, only: ArModel
   implicit none
contains
   subroutine fugacity_tp(self, &
         n, T, P, V, root_type, lnfug, dlnPhidP, dlnphidT, dlnPhidn &
      )
      use iso_fortran_env, only: error_unit
      class(ArModel), intent(in) :: self
      real(pr), intent(in) :: n(:) !! Mixture mole numbers
      integer, intent(in) :: root_type !! Type of root desired (-1 vapor, 1 liquid, 0 lower Gr)
      real(pr), intent(in) :: t    !! Temperature [K]
      real(pr), intent(in) :: p    !! Pressure [bar]

      real(pr), intent(out) :: lnfug(size(n)) !! \(\ln(\phi*p)\) vector
      real(pr), optional, intent(out) :: v !! Volume [L]
      real(pr), optional, intent(out) :: dlnphidt(size(n)) !! ln(phi) Temp derivative
      real(pr), optional, intent(out) :: dlnphidp(size(n)) !! ln(phi) Presssure derivative
      real(pr), optional, intent(out) :: dlnphidn(size(n), size(n)) !! ln(phi) compositional derivative

      real(pr) :: v_in, p_in

      call VCALC(self, root_type, size(n), n, T, P, V_in)
      call fugacity_vt(self, n, v_in, T, P_in, lnfug, dlnphidp, dlnphidt, dlnphidn)
      
      if(present(v)) v = v_in
      if(abs(p-p_in) > 1e-5) write(error_unit, *) &
       "WARN: Possible wrong root calculating fugacity!" , p, p_in
   end subroutine fugacity_tp

   subroutine fugacity_vt(self, &
         n, V, T, P, lnfug, dlnPhidP, dlnphidT, dlnPhidn &
      )
      class(ArModel) :: self
      real(pr), intent(in) :: n(:) !! Mixture mole numbers
      real(pr), intent(in) :: v !! Volume [L]
      real(pr), intent(in) :: t !! Temperature [K]

      real(pr), optional, intent(out) :: p    !! Pressure [bar]
      real(pr), optional, intent(out) :: lnfug(size(n)) !! \(\ln(\phi*p)\) vector
      real(pr), optional, intent(out) :: dlnphidt(size(n)) !! ln(phi) Temp derivative
      real(pr), optional, intent(out) :: dlnphidp(size(n)) !! ln(phi) Presssure derivative
      real(pr), optional, intent(out) :: dlnphidn(size(n), size(n)) !! ln(phi) compositional derivative

      real(pr) :: Ar, ArTV, ArV, ArV2
      real(pr), dimension(size(n)) :: Arn, ArVn, ArTn
      real(pr) :: Arn2(size(n), size(n))

      real(pr) :: dPdV, dPdT, dPdN(size(n))

      real(pr) :: RT, Z

      real(pr) :: totn
      integer :: nc, i, j


      TOTN = sum(n)
      nc = size(n)

      RT = R*T
      Z = V/(TOTN*RT) ! this is Z/P

      !if (present(dlnphidn)) then
         call self%residual_helmholtz(&
            n, V, T, Ar=Ar, ArV=ArV, ArV2=ArV2, ArTV=ArTV, &
            Arn=Arn, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2 &
         )
      ! else
      !    call self%residual_helmholtz(&
      !       n, V, T, Ar=Ar, ArV=ArV, ArV2=ArV2, ArTV=ArTV, &
      !       Arn=Arn, ArVn=ArVn, ArTn=ArTn &
      !    )
      ! end if

      P = TOTN*RT/V - ArV
      dPdV = -ArV2 - RT*TOTN/V**2
      dPdT = -ArTV + TOTN*R/V
      dPdN(:) = RT/V - ArVn(:)

      lnfug(:) = Arn(:)/RT - log(Z)

      if (present(dlnphidp)) dlnphidp(:) = -dPdN(:)/dPdV/RT - 1.D0/P
      if (present(dlnphidt)) then
         dlnphidt(:) = (ArTn(:) - Arn(:)/T)/RT + dPdN(:)*dPdT/dPdV/RT + 1.D0/T
      end if

      if (present(dlnphidn)) then
         do i = 1, nc
            do j = i, nc
               dlnphidn(i, j) = 1.D0/TOTN + (Arn2(i, j) + dPdN(i)*dPdN(j)/dPdV)/RT
               dlnphidn(j, i) = dlnphidn(i, j)
            end do
         end do
      end if
   end subroutine

   subroutine PUREFUG_CALC(self, nc, icomp, T, P, V, fug)
      !! Fugacity of a pure component
      class(ArModel), intent(in) :: self
      integer, intent(in)  :: nc
      integer,  intent(in) :: icomp
      real(pr), intent(in) :: T, P, V
      real(pr), intent(out) :: fug

      real(pr) :: n(nc), Ar, Arn(nc)
      real(pr) :: RT, Z, lnfug

      n = 0.0
      n(icomp) = 1.0

      RT = R*T
      Z = P*V/RT

      call self%residual_helmholtz(n, V, T, Ar, Arn=Arn)
      lnfug = -log(Z) + Arn(icomp)/RT
      fug = exp(fug)
   end subroutine purefug_calc

   recursive subroutine VCALC(self, ITYP, nc, rn, T, P, V, max_iters)
      !! ROUTINE FOR CALCULATION OF VOLUME, GIVEN PRESSURE
      use iso_fortran_env, only: error_unit
      use stdlib_optval, only: optval
      class(ArModel), intent(in) :: self
      integer, intent(in) :: ITYP !! Type of root [-1: vapor, 1:liquid, 0:lower Gibbs energy phase]
      integer, intent(in) :: nc  !! Number of components
      real(pr), intent(in) ::  rn(nc) !! Mixture moles
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(out) :: V !! Volume [L]
      integer, optional, intent(in) :: max_iters !! Maxiumum number of iterations, defaults to 100

      real(pr) ::  Ar, ArV, ArV2
      logical :: FIRST_RUN

      real(pr) :: totn
      real(pr) :: B, CPV, S3R
      real(pr) :: ZETMIN, ZETA, ZETMAX
      real(pr) :: del, pcalc, der, AT, AVAP, VVAP

      integer :: iter, maximum_iterations

      maximum_iterations = optval(max_iters, 100)

      FIRST_RUN = .true.
      TOTN = sum(rn)
      CPV = self%get_v0(rn, p, t)
      B = CPV
      S3R = 1.D0/CPV
      ITER = 0

      ZETMIN = 0.D0
      ZETMAX = 1.D0 - 0.01*T/(10000*B)  ! improvement for cases with heavy components
      if (ITYP .gt. 0) then
         ZETA = .5D0
      else
         ! IDEAL GAS ESTIMATE
         ZETA = min(.5D0, CPV*P/(TOTN*R*T))
      end if
  100 continue
      DEL = 1
      pcalc = 2*p

      do while(abs(DEL) > 1d-10 .and. iter < maximum_iterations)
         V = CPV/ZETA
         ITER = ITER + 1
         call self%residual_helmholtz(&
            rn, V, T, Ar=Ar, ArV=ArV, ArV2=ArV2 &
         )
         PCALC = TOTN*R*T/V - ArV

         if (PCALC .gt. P) then
            ZETMAX = ZETA
         else
            ZETMIN = ZETA
         end if

         ! AT is something close to Gr(P,T)
         AT = (Ar + V*P)/(T*R) - TOTN*log(V)

         DER = (ArV2*V**2 + TOTN*R*T)*S3R  ! this is dPdrho/B
         DEL = -(PCALC - P)/DER
         ZETA = ZETA + max(min(DEL, 0.1D0), -.1D0)

         if (ZETA .gt. ZETMAX .or. ZETA .lt. ZETMIN) &
            ZETA = .5D0*(ZETMAX + ZETMIN)
      end do

      if (iter >= maximum_iterations) write(error_unit, *) &
         "WARN: Volume solver exceeded maximum number of iterations"

      if (ITYP .eq. 0) then
         ! FIRST RUN WAS VAPOUR; RERUN FOR LIQUID
         if (FIRST_RUN) then
            VVAP = V
            AVAP = AT
            FIRST_RUN = .false.
            ZETA = 0.5D0
            ZETMAX = 1.D0 - 0.01*T/500
            goto 100
         else
            if (AT .gt. AVAP) V = VVAP
         end if
      end if
   end subroutine vcalc
   ! ==========================================================================
end module