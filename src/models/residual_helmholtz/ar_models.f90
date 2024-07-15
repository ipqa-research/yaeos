module yaeos__models_ar
  
   !! # Module that defines the basics of a residual Helmholtz energy.
   !!
   !! All the residual properties that are calculated in this library are
   !! based on residual Helmholtz Equations of State. Following the book by
   !! Michelsen and Mollerup.
   !!
   !! In this library up to second derivatives of residual Helmholtz energy
   !! are used. Because they're the fundamentals for phase equilibria
   !! calculation.
   !!
   !! @note 
   !! Later on, third derivative with respect to volume will be included
   !! since it's importance on calculation of critical points.
   !! @endnote
   !!
   !! # Properties
   !!
   !! ## Available properties:
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
   use yaeos__constants, only: pr, R
   use yaeos__models_base, only: BaseModel
   implicit none

   type, abstract, extends(BaseModel) :: ArModel
      !! Abstract residual Helmholtz model.
      !!
      !! This derived type defines the basics needed for the calculation
      !! of residual properties.
      !! The basics of a residual Helmholtz model is a routine that calculates
      !! all the needed derivatives of \(Ar\) `residual_helmholtz` and
      !! a volume initializer function, that is used to initialize a Newton
      !! solver of volume when specifying pressure.
      character(len=:), allocatable :: name !! Name of the model
   contains
      procedure(abs_residual_helmholtz), deferred :: residual_helmholtz
      !! Method to calculate residual helmholtz energy and derivatives.
      procedure(abs_volume_initializer), deferred :: get_v0
      !! Volume initializer
      procedure :: lnphi_vt => fugacity_vt
      procedure :: lnphi_tp => fugacity_tp
      procedure :: pressure
      procedure :: volume
      procedure :: enthalpy_residual_vt
      procedure :: gibbs_residual_vt
      procedure :: entropy_residual_vt
      procedure :: Cv_residual_vt
      procedure :: Cp_residual_vt
   end type ArModel

   abstract interface
      subroutine abs_residual_helmholtz(&
         self, n, v, t, Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2 &
         )
         !! Residual Helmholtz model generic interface.
         !!
         !! This interface represents how an Ar model should be implemented.
         !! By our standard, a Resiudal Helmholtz model takes as input:
         !!
         !! - The mixture's number of moles vector.
         !! - Volume, by default in liters.
         !! - Temperature, by default in Kelvin.
         !!
         !! All the output arguments are optional. While this keeps a long
         !! signature for the implementation, this is done this way to take
         !! advantage of any inner optimizations to calculate derivatives
         !! inside the procedure.
         !!
         !! Once the model is implemented, the signature can be short like
         !! `model%residual_helmholtz(n, v, t, ArT2=dArdT2)`
         import ArModel, pr
         class(ArModel), intent(in) :: self !! ArModel
         real(pr), intent(in) :: n(:) !! Moles vector
         real(pr), intent(in) :: v !! Volume [L]
         real(pr), intent(in) :: t !! Temperature [K]
         real(pr), optional, intent(out) :: Ar !! Residual Helmoltz energy
         real(pr), optional, intent(out) :: ArV !! \(\frac{dAr}{dV}\)
         real(pr), optional, intent(out) :: ArT !! \(\frac{dAr}{dT}\)
         real(pr), optional, intent(out) :: ArT2 !! \(\frac{d^2Ar}{dT^2}\)
         real(pr), optional, intent(out) :: ArTV !! \(\frac{d^2Ar}{dTV}\)
         real(pr), optional, intent(out) :: ArV2 !! \(\frac{d^2Ar}{dV^2}\)
         real(pr), optional, intent(out) :: Arn(size(n)) !! \(\frac{dAr}{dn_i}\)
         real(pr), optional, intent(out) :: ArVn(size(n)) !! \(\frac{d^2Ar}{dVn_i}\)
         real(pr), optional, intent(out) :: ArTn(size(n)) !! \(\frac{d^2Ar}{dTn_i}\)
         real(pr), optional, intent(out) :: Arn2(size(n), size(n))!! \(\frac{d^2Ar}{dn_{ij}}\)
      end subroutine abs_residual_helmholtz

      function abs_volume_initializer(self, n, p, t)
         !! Initialization of volume.
         import ArModel, pr
         class(ArModel), intent(in) :: self !! Ar Model
         real(pr), intent(in) :: n(:) !! Moles vector
         real(pr), intent(in) :: p !! Pressure [bar]
         real(pr), intent(in) :: t !! Temperature [K]
         real(pr) :: abs_volume_initializer !! Initial volume [L]
      end function abs_volume_initializer
   end interface

contains

   subroutine volume(eos, n, P, T, V, root_type)
      !! Generic volume solver
      use yaeos__constants, only: pr, R
      use yaeos__math, only: newton
      class(ArModel), intent(in) :: eos
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: P, T
      real(pr), intent(out) :: V
      character(len=*), intent(in) :: root_type

      integer :: max_iters=30
      real(pr) :: tol=1e-7

      real(pr) :: totnRT, GrL, GrV, Gr
      real(pr) :: Vliq, Vvap

      GrL = HUGE(GrL)
      GrV = HUGE(GrV)

      totnRT = sum(n) * R * T
      select case(root_type)
       case("liquid")
         Vliq = eos%get_v0(n, P, T)! *1.001_pr
         call newton(foo, Vliq, tol=tol, max_iters=max_iters)
         GrL = Gr
       case("vapor")
         Vvap = R * T / P
         call newton(foo, Vvap, tol=tol, max_iters=max_iters)
         GrV = Gr
       case("stable")
         Vliq = eos%get_v0(n, P, T)*1.00001_pr
         call newton(foo, Vliq, tol=tol, max_iters=max_iters)
         GrL = Gr

         Vvap = R * T / P
         call newton(foo, Vvap, tol=tol, max_iters=max_iters)
         GrV = Gr
      end select

      if (GrL < GrV) then
         V = Vliq
      else
         V = Vvap
      end if

   contains
      subroutine foo(x, f, df)
         real(pr), intent(in) :: x
         real(pr), intent(out) :: f, df
         real(pr) :: Ar, ArV, ArV2, Pcalc, dPcalcdV, Vin
         Vin = x
         call eos%residual_helmholtz(n, Vin, T, Ar=Ar, ArV=ArV, ArV2=ArV2)
         Pcalc = totnRT / Vin - ArV
         dPcalcdV = -totnRT / Vin**2 - ArV2
         f = Pcalc - p
         df = dPcalcdV
         Gr = Ar + P * Vin - totnRT - totnRT * log(P*Vin/(R*T))
      end subroutine foo
   end subroutine volume

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
      logical :: dn

      TOTN = sum(n)
      nc = size(n)

      if (present(dPdN)) then
         call eos%residual_helmholtz(&
            n, v, t, Ar=Ar, ArV=ArV, ArV2=ArV2, ArTV=ArTV, ArVn=ArVn &
            )
      else
         call eos%residual_helmholtz(&
            n, v, t, Ar=Ar, ArV=ArV, ArV2=ArV2, ArTV=ArTV &
            )
      end if

      P = TOTN*R*T/V - ArV
      if (present(dPdV)) dPdV = -ArV2 - R*t*TOTN/V**2
      if (present(dPdT)) dPdT = -ArTV + TOTN*R/V
      if (present(dPdN)) dPdN(:) = R*T/V - ArVn(:)
   end subroutine pressure

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
      if(abs(P_in - p) > 1e-2) then
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
      real(pr) :: P_in

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

      P_in = TOTN*RT/V - ArV
      if (present(P)) P = P_in

      dPdV_in = -ArV2 - RT*TOTN/V**2
      dPdT_in = -ArTV + TOTN*R/V
      dPdN_in = RT/V - ArVn

      if (present(dlnphidp)) then
         dlnphidp(:) = -dPdN_in(:)/dPdV_in/RT - 1._pr/P_in
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
   end subroutine fugacity_vt

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
         n, v, t, Ar=Ar, ArV=ArV, ArT=ArT, ArTV=ArTV, ArV2=ArV2, ArT2=ArT2, Arn=Arn, ArVn=ArVn, ArTn=ArTn &
         )

      Hr = Ar - t*ArT - v*ArV

      if (present(HrT)) HrT = - t*ArT2 - v*ArTV
      if (present(HrV)) HrV = - t*ArTV - v*ArV2
      if (present(HrN)) HrN(:) = Arn(:) - t*ArTn(:) - v*ArVn(:)
   end subroutine enthalpy_residual_vt

   subroutine gibbs_residual_vt(eos, n, v, t, Gr, GrT, GrV, Grn)
      !! Calculate residual Gibbs energy given volume and temperature.
      use yaeos__constants, only: R
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
   end subroutine gibbs_residual_vt

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
   end subroutine entropy_residual_vt

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
   end subroutine Cv_residual_vt

   subroutine Cp_residual_vt(eos, n, v, t, Cp)
      !! Calculate residual heat capacity pressure constant given v and t.
      use yaeos__constants, only: R
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
   end subroutine Cp_residual_vt
end module yaeos__models_ar
