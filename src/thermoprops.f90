module yaeos_thermoprops
   !! Residual thermodyamic properties using residual Helmholtz model.
   !! 
   !! Available properties:
   !!
   !! - pressure(n, V, T)
   !! - fugacity(n, V, T)
   !! - fugacity(n, P, T, root=[vapor, liquid, stable])
   !! - volume
   !!
   !! Calculate thermodynamic properties using Helmholtz energy as a basis.
   !! All the routines in this module work with the logic:
   !! 
   !! ```fortran
   !! call foo(x, V, T, [dfoodv, dfoodt, ...])
   !! ```
   !! Where the user can call the routine of the desired property. And include 
   !! as optional values the desired derivatives of said properties.
   use yaeos_constants, only: R, pr
   use yaeos_models_ar, only: ArModel
   implicit none
contains
   subroutine pressure(eos, n, v, t, p, dpdv, dpdt, dpdn)
      !! Pressure calculation.
      !!
      !! Calculate pressure using residual helmholtz models.
      !!
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector 
      real(pr), intent(in) :: t !! Temperature [K]
      real(pr), intent(in) :: v !! Volume [L]
      real(pr), intent(out) :: p !! Pressure [bar]
      real(pr), optional, intent(out) :: dpdv !! \(\frac{dP}}{dV}\)
      real(pr), optional, intent(out) :: dpdt !! \(\frac{dP}}{dT}\)
      real(pr), optional, intent(out) :: dpdn(:) !! \(\frac{dP}}{dn_i}\)

      real(pr) :: totn

      real(pr) :: Ar, ArV, ArV2, ArTV, ArVn(size(n))
      integer :: nc
      
      TOTN = sum(n)
      nc = size(n)

      call eos%residual_helmholtz(&
         n, v, t, Ar=Ar, ArV=ArV, ArV2=ArV2, ArTV=ArTV, ArVn=ArVn &
      )
      P = TOTN*R*T/V - ArV
      if (present(dPdV)) dPdV = -ArV2 - R*t*TOTN/V**2
      if (present(dPdT)) dPdT = -ArTV + TOTN*R/V
      if (present(dPdN)) dPdN(:) = R*T/V - ArVn(:)
   end subroutine

   subroutine fugacity_tp(eos, &
         n, T, P, V, root_type, lnphip, dlnPhidP, dlnphidT, dlnPhidn &
      )
      !! Calculate logarithm of fugacity, given pressure and temperature.
      !! 
      !! This routine will obtain the desired volume root at the specified
      !! pressure and calculate fugacity at that point.
      !!
      !! @note
      !! While the natural output variable is \(ln f_i\). The calculated
      !! derivatives will be the derivatives of the fugacity coefficient
      !! \(ln \phi_i\)
      !! @endnote
      !!
      use iso_fortran_env, only: error_unit
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Mixture mole numbers
      character(len=*), intent(in) :: root_type !! Type of root desired ["liquid", "vapor", "stable"]
      real(pr), intent(in) :: t    !! Temperature [K]
      real(pr), intent(in) :: p    !! Pressure [bar]

      real(pr), optional, intent(out) :: lnphip(size(n)) !! \(\ln(f_i)\) vector
      real(pr), optional, intent(out) :: v !! Volume [L]
      real(pr), optional, intent(out) :: dlnphidt(size(n)) !! ln(phi) Temp derivative
      real(pr), optional, intent(out) :: dlnphidp(size(n)) !! ln(phi) Presssure derivative
      real(pr), optional, intent(out) :: dlnphidn(size(n), size(n)) !! ln(phi) compositional derivative

      real(pr) :: v_in, p_in

      call volume(eos, n, P, T, V_in, root_type)
      call fugacity_vt(eos, n, v_in, T, P_in, lnphip, dlnphidp, dlnphidt, dlnphidn)
      if(present(v)) v = v_in
      if(abs(P_in - p)/p > 1e-5) then
         write(error_unit, *) "WARN: possible bad root solving: ", p_in, p
      end if
   end subroutine fugacity_tp

   subroutine fugacity_vt(eos, &
         n, V, T, P, lnphip, dlnPhidP, dlnphidT, dlnPhidn, dPdV, dPdT, dPdN &
      )
      !! Calculate fugacity coefficent given volume and temperature.
      !!
      !!@note
      !!While the natural output variable is \(ln \phi_i P\). The calculated
      !!derivatives will be the derivatives of the fugacity coefficient
      !!\(ln \phi_i\)
      !!@endnote
      !!
      class(ArModel) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Mixture mole numbers
      real(pr), intent(in) :: v !! Volume [L]
      real(pr), intent(in) :: t !! Temperature [K]

      real(pr), optional, intent(out) :: p !! Pressure [bar]
      real(pr), optional, intent(out) :: lnphip(size(n)) !! \(\ln(\phi_i*P)\) vector
      real(pr), optional, intent(out) :: dlnphidt(size(n)) !! \(ln(phi_i)\) Temp derivative
      real(pr), optional, intent(out) :: dlnphidp(size(n)) !! \(ln(phi_i)\) Presssure derivative
      real(pr), optional, intent(out) :: dlnphidn(size(n), size(n)) !! \(ln(phi_i)\) compositional derivative
      real(pr), optional, intent(out) :: dPdV !! \(\frac{dP}{dV}\)
      real(pr), optional, intent(out) :: dPdT !! \(\frac{dP}{dT}\)
      real(pr), optional, intent(out) :: dPdN(:) !! \(\frac{dP}{dn_i}\)

      real(pr) :: Ar, ArTV, ArV, ArV2
      real(pr), dimension(size(n)) :: Arn, ArVn, ArTn
      real(pr) :: Arn2(size(n), size(n))

      real(pr) :: dPdV_in, dPdT_in, dPdN_in(size(n))

      real(pr) :: RT, Z

      real(pr) :: totn
      integer :: nc, i, j


      totn = sum(n)
      nc = size(n)

      RT = R*T
      Z = V/(TOTN*RT) ! this is Z/P

      if (present(lnphip) .and. .not. (&
              present(dlnphidn) &
         .or. present(dlnphidp) &
         .or. present(dlnphidt) &
         .or. present(p) &
      )) then
         call eos%residual_helmholtz(n, v, t, Arn=Arn)
         lnphip(:) = Arn(:)/RT - log(Z)
         return
      else if (present(dlnphidn)) then
         call eos%residual_helmholtz(&
            n, V, T, Ar=Ar, ArV=ArV, ArV2=ArV2, ArTV=ArTV, &
            Arn=Arn, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2 &
         )
      else
         call eos%residual_helmholtz(&
            n, V, T, Ar=Ar, ArV=ArV, ArV2=ArV2, ArTV=ArTV, &
            Arn=Arn, ArVn=ArVn, ArTn=ArTn &
         )
      end if
      
      lnphip(:) = Arn(:)/RT - log(Z)

      if (present(P)) P = TOTN*RT/V - ArV

      dPdV_in = -ArV2 - RT*TOTN/V**2
      dPdT_in = -ArTV + TOTN*R/V
      dPdN_in = RT/V - ArVn

      if (present(dlnphidp)) then
         dlnphidp(:) = -dPdN_in(:)/dPdV_in/RT - 1._pr/P
      end if
      if (present(dlnphidt)) then
         dlnphidt(:) = (ArTn(:) - Arn(:)/T)/RT + dPdN_in(:)*dPdT_in/dPdV_in/RT + 1._pr/T
      end if

      if (present(dlnphidn)) then
         do i = 1, nc
            do j = i, nc
               dlnphidn(i, j) = 1._pr/TOTN + (Arn2(i, j) + dPdN_in(i)*dPdN_in(j)/dPdV_in)/RT
               dlnphidn(j, i) = dlnphidn(i, j)
            end do
         end do
      end if

      if (present(dPdV)) dPdV = dPdV_in
      if (present(dPdT)) dPdT = dPdT_in
      if (present(dPdN)) dPdN = dPdN_in
   end subroutine

   subroutine PUREFUG_CALC(eos, nc, icomp, T, P, V, fug)
      !! Fugacity of a pure component.
      class(ArModel), intent(in) :: eos !! model
      integer, intent(in)  :: nc !! Number of components
      integer,  intent(in) :: icomp !! Component index
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: V !! Volume [L]
      real(pr), intent(out) :: fug !! Fugacity of component `icomp`

      real(pr) :: n(nc), Ar, Arn(nc)
      real(pr) :: RT, Z, lnphi

      n = 0.0
      n(icomp) = 1.0

      RT = R*T
      Z = P*V/RT

      call eos%residual_helmholtz(n, V, T, Ar, Arn=Arn)
      lnphi = -log(Z) + Arn(icomp)/RT
      fug = exp(lnphi) * P
   end subroutine purefug_calc

   subroutine volume(eos, n, P, T, V, root_type, max_iters)
      !! Volume solver at a given pressure.
      !!
      !! Obtain the volume using the method described by Michelsen and Møllerup.
      !! While \(P(V, T)\) can be obtained with a simple Newton method, a better
      !! approach is solving \(P(B/V, T)\) where \(B\) is the EoS covolume.
      !! This method is easier to solve because:
      !! \[
      !!    V(P, T) \in [0, \infty)
      !! \]
      !! and
      !! \[
      !!    \frac{B}{V}(P, T) \in [0, 1]
      !! \]
      !! 
      !! At chapter 3 page 94 of Michelsen and Møllerup's book a more complete
      !! explanation can be seen
      use iso_fortran_env, only: error_unit
      use stdlib_optval, only: optval
      class(ArModel), intent(in) :: eos
      real(pr), intent(in) ::  n(:) !! Mixture moles
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(out) :: V !! Volume [L]
      character(len=*), optional, intent(in) :: root_type !! Type of root ["vapor" | "liquid" | "stable"]
      integer, optional, intent(in) :: max_iters !! Maxiumum number of iterations, defaults to 100

      character(len=10) :: root

      real(pr) ::  Ar, ArV, ArV2

      real(pr) :: totn
      real(pr) :: B !! Covolume
      real(pr) :: ZETMIN, ZETA, ZETMAX
      real(pr) :: pcalc, AT, AVAP, VVAP

      integer :: iter, maximum_iterations

      maximum_iterations = optval(max_iters, 100)
      root = optval(root_type, "stable")

      TOTN = sum(n)
      B = eos%get_v0(n, p, t)
      ITER = 0

      ! Limits
      ZETMIN = 0._pr
      ZETMAX = 1._pr

      select case(root_type)
      case("liquid")
         ZETA = 0.5_pr
         call solve_point(P, V, Pcalc, AT, iter)
      case("vapor","stable")
         ZETA = min(0.5_pr, B*P/(TOTN*R*T))
         call solve_point(P, V, Pcalc, AT, iter)
         
         if (root_type == "stable") then
            ! Run first for vapor and then for liquid
            VVAP = V
            AVAP = AT
            ZETA = 0.5_pr
            ZETMAX = 1._pr
            call solve_point(P, V, Pcalc, AT, iter)
            if (AT .gt. AVAP) V = VVAP
         end if
      case default
         write(error_unit, *) "ERROR [VCALC]: Wrong specification"
         error stop 1
      end select
   contains
      subroutine solve_point(P, V, Pcalc, AT, iter)
         real(pr), intent(in) :: P !! Objective pressure [bar]
         real(pr), intent(out) :: V !! Obtained volume [L]
         real(pr), intent(out) :: Pcalc !! Calculated pressure at V [bar]
         real(pr), intent(out) :: AT !!
         integer, intent(out) :: iter 

         real(pr) :: del, der

         iter = 0
         DEL = 1
         pcalc = 2*p
         do while(abs(DEL) > 1.e-10_pr .and. iter < maximum_iterations)
            V = B/ZETA
            iter = iter + 1
            call eos%residual_helmholtz(n, V, T, Ar=Ar, ArV=ArV, ArV2=ArV2)
            
            Pcalc = TOTN*R*T/V - ArV

            if (Pcalc .gt. P) then
               ZETMAX = ZETA
            else
               ZETMIN = ZETA
            end if

            ! AT is something close to Gr(P,T)
            AT = (Ar + V*P)/(T*R) - TOTN*log(V)

            DER = (ArV2*V**2 + TOTN*R*T)/B  ! this is dPdrho/B
            DEL = -(Pcalc - P)/DER
            ZETA = ZETA + max(min(DEL, 0.1_pr), -.1_pr)

            if (ZETA .gt. ZETMAX .or. ZETA .lt. ZETMIN) then 
               ZETA = 0.5_pr*(ZETMAX + ZETMIN)
            end if
         end do
         
         if (iter >= maximum_iterations) write(error_unit, *) &
            "WARN: Volume solver exceeded maximum number of iterations"
      end subroutine
   end subroutine volume
   
   ! ==========================================================================
   ! Residual Enthalpy
   ! --------------------------------------------------------------------------
   subroutine enthalpy_residual_vt(eos, n, v, t, Hr, HrT, HrV, Hrn)
      !! Calculate residual enthalpy given volume and temperature.
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: t !! Temperature [K]
      real(pr), intent(in) :: v !! Volume [L]
      real(pr), intent(out) :: Hr !! Residual enthalpy [bar L / mol]
      real(pr), optional, intent(out) :: HrT !! \(\frac{dH^r}}{dT}\)
      real(pr), optional, intent(out) :: HrV !! \(\frac{dH^r}}{dV}\)
      real(pr), optional, intent(out) :: Hrn(size(n)) !! \(\frac{dH^r}}{dn}\)

      real(pr) :: Ar, ArV, ArT, Arn(size(n))
      real(pr) :: ArV2, ArT2, ArTV, ArVn(size(n)), ArTn(size(n))

      call eos%residual_helmholtz(&
         n, v, t, Ar=Ar, ArV=ArV, ArT=ArT, Arn=Arn, &
         ArTV=ArTV, ArV2=ArV2, ArT2=ArT2, ArVn=ArVn, ArTn=ArTn &
      )

      Hr = Ar - t*ArT - v*ArV

      if (present(HrT)) HrT = - t*ArT2 - v*ArTV
      if (present(HrV)) HrV = - t*ArTV - v*ArV2
      if (present(HrN)) HrN(:) = Arn(:) - t*ArTn(:) - v*ArVn(:)
   end subroutine

   ! ==========================================================================
   ! Residual Gibbs energy
   ! --------------------------------------------------------------------------
   subroutine gibbs_residual_vt(eos, n, v, t, Gr, GrT, GrV, Grn)
      !! Calculate residual Gibbs energy given volume and temperature.
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: t !! Temperature [K]
      real(pr), intent(in) :: v !! Volume [L]
      real(pr), intent(out) :: Gr !! Gibbs energy [bar L / mol]
      real(pr), optional, intent(out) :: GrT !! \(\frac{dG^r}}{dT}\)
      real(pr), optional, intent(out) :: GrV !! \(\frac{dG^r}}{dV}\)
      real(pr), optional, intent(out) :: Grn(size(n)) !! \(\frac{dG^r}}{dn}\)

      real(pr) :: Ar, ArV, ArT, Arn(size(n))
      real(pr) :: p, dpdv, dpdt, dpdn(size(n)), z, ntot

      ntot = sum(n)
      call pressure(eos, n, v, t, p, dpdv=dpdv, dpdt=dpdt, dpdn=dpdn)
      z = p*v/(ntot*R*t)

      call eos%residual_helmholtz(n, v, t, Ar=Ar, ArV=ArV, ArT=ArT, Arn=Arn)

      Gr = Ar + p*v - ntot*R*t

      if (present(GrT)) GrT = ArT + v*dpdt - ntot*R
      if (present(GrV)) GrV = ArV + v*dpdv + p
      if (present(GrN)) GrN(:) = Arn(:) + v*dpdn(:) - R*t
   end subroutine

   ! ==========================================================================
   ! Residual entropy
   ! --------------------------------------------------------------------------
   subroutine entropy_residual_vt(eos, n, v, t, Sr, SrT, SrV, Srn)
      !! Calculate residual entropy given volume and temperature.
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: t !! Temperature [K]
      real(pr), intent(in) :: v !! Volume [L]
      real(pr), intent(out) :: Sr !! Entropy [bar L / K / mol]
      real(pr), optional, intent(out) :: SrT !! \(\frac{dS^r}}{dT}\)
      real(pr), optional, intent(out) :: SrV !! \(\frac{dS^r}}{dV}\)
      real(pr), optional, intent(out) :: Srn(size(n)) !! \(\frac{dS^r}}{dn}\)

      real(pr) :: Ar, ArT, ArT2, ArTV, ArTn(size(n))

      call eos%residual_helmholtz(&
         n, v, t, Ar=Ar, ArT=ArT, ArTV=ArTV, ArT2=ArT2, ArTn=ArTn &
      )

      Sr = - ArT

      if (present(SrT)) SrT = - ArT2
      if (present(SrV)) SrV = - ArTV
      if (present(SrN)) SrN = - ArTn
   end subroutine

   ! ==========================================================================
   ! Residual Cv
   ! --------------------------------------------------------------------------
   subroutine Cv_residual_vt(eos, n, v, t, Cv)
      !! Calculate residual heat capacity volume constant given v and t.
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: t !! Temperature [K]
      real(pr), intent(in) :: v !! Volume [L]
      real(pr), intent(out) :: Cv !! heat capacity v constant [bar L / K / mol]

      real(pr) :: Ar, ArT2

      call eos%residual_helmholtz(n, v, t, Ar=Ar, ArT2=ArT2)

      Cv = -t*ArT2
   end subroutine

   ! ==========================================================================
   ! Residual Cp
   ! --------------------------------------------------------------------------
   subroutine Cp_residual_vt(eos, n, v, t, Cp)
      !! Calculate residual heat capacity pressure constant given v and t.
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: t !! Temperature [K]
      real(pr), intent(in) :: v !! Volume [L]
      real(pr), intent(out) :: Cp !! heat capacity p constant [bar L / K / mol]

      real(pr) :: Ar, ArT2, Cv, p, dpdt, dpdv, ntot

      ntot = sum(n)

      call eos%residual_helmholtz(n, v, t, Ar=Ar, ArT2=ArT2)

      call Cv_residual_vt(eos, n, v, t, Cv)

      call pressure(eos, n, v, t, p, dpdv=dpdv, dpdt=dpdt)

      Cp = Cv - t*dpdt**2/dpdv - ntot*R
   end subroutine
end module
