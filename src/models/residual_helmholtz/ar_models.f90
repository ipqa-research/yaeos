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
   !! call foo(x, V, T, [dfoodv, dfoodT, ...])
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
      procedure(abs_volume_initializer), deferred :: get_v0
      procedure :: lnphi_vt => fugacity_vt
      procedure :: lnphi_pt => fugacity_pt
      procedure :: pressure
      procedure :: volume
      procedure :: enthalpy_residual_vt
      procedure :: gibbs_residual_vt
      procedure :: entropy_residual_vt
      procedure :: Cv_residual_vt
      procedure :: Cp_residual_vt
   end type ArModel

   interface size
      module procedure :: size_ar_model
   end interface size

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
         !! Function that provides an initializer value for the liquid-root
         !! of newton solver of volume. In the case the model will use the
         !! `volume_michelsen` routine this value should provide the co-volume
         !! of the model.
         import ArModel, pr
         class(ArModel), intent(in) :: self !! Ar Model
         real(pr), intent(in) :: n(:) !! Moles vector
         real(pr), intent(in) :: p !! Pressure [bar]
         real(pr), intent(in) :: t !! Temperature [K]
         real(pr) :: abs_volume_initializer !! Initial volume [L]
      end function abs_volume_initializer
   end interface

contains

   integer pure function size_ar_model(eos)
      !! Get the size of the model.
      class(ArModel), intent(in) :: eos
      size_ar_model = size(eos%components%pc)
   end function size_ar_model

   subroutine volume(eos, n, P, T, V, root_type)
      !! # Volume solver routine for residual Helmholtz models.
      !! Solves volume roots using newton method. Given pressure and temperature.
      !!
      !! # Description
      !! This subroutine solves the volume using a newton method. The variable
      !! `root_type`
      !!
      !! # Examples
      !!
      !! ```fortran
      !! class(ArModel) :: eos
      !! call eos%volume(n, P, T, V, root_type="liquid")
      !! call eos%volume(n, P, T, V, root_type="vapor")
      !! call eos%volume(n, P, T, V, root_type="stable")
      !! ```
      use yaeos__constants, only: pr, R
      use yaeos__math, only: newton
      class(ArModel), intent(in) :: eos
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(out) :: V !! Volume [L]
      character(len=*), intent(in) :: root_type
      !! Desired root-type to solve. Options are:
      !! `["liquid", "vapor", "stable"]`

      integer :: max_iters=30
      real(pr) :: tol=1e-7

      real(pr) :: totnRT, GrL, GrV, Gr
      real(pr) :: Vliq, Vvap

      GrL = HUGE(GrL)
      GrV = HUGE(GrV)

      totnRT = sum(n) * R * T
      select case(root_type)
       case("liquid")
         Vliq = eos%get_v0(n, P, T) * 1.001_pr
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

   subroutine pressure(eos, n, V, T, P, dPdV, dPdT, dPdn)
      !! Pressure calculation.
      !!
      !! Calculate pressure using residual helmholtz models.
      !!
      !! # Examples
      !! ```fortran
      !! class(ArModel), allocatable :: eos
      !! real(pr) :: n(2), t, v, p, dPdV, dPdT, dPdn(2)
      !! eos = PengRobinson(Tc, Pc, w)
      !! n = [1.0_pr, 1.0_pr]
      !! t = 300.0_pr
      !! v = 1.0_pr
      !! call eos%pressure(n, V, T, P, dPdV=dPdV, dPdT=dPdT, dPdn=dPdn)
      !! ```
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: t !! Temperature [K]
      real(pr), intent(in) :: v !! Volume [L]
      real(pr), intent(out) :: p !! Pressure [bar]
      real(pr), optional, intent(out) :: dPdV !! \(\frac{dP}}{dV}\)
      real(pr), optional, intent(out) :: dPdT !! \(\frac{dP}}{dT}\)
      real(pr), optional, intent(out) :: dPdn(:) !! \(\frac{dP}}{dn_i}\)

      real(pr) :: totn

      real(pr) :: Ar, ArV, ArV2, ArTV, ArVn(size(eos))
      integer :: nc
      logical :: dn

      totn = sum(n)
      nc = size(n)

      if (present(dPdn)) then
         call eos%residual_helmholtz(&
            n, v, t, Ar=Ar, ArV=ArV, ArV2=ArV2, ArTV=ArTV, ArVn=ArVn &
            )
      else
         call eos%residual_helmholtz(&
            n, v, t, Ar=Ar, ArV=ArV, ArV2=ArV2, ArTV=ArTV &
            )
      end if

      P = totn * R * T/V - ArV
      if (present(dPdV)) dPdV = -ArV2 - R*T*totn/V**2
      if (present(dPdT)) dPdT = -ArTV + totn*R/V
      if (present(dPdn)) dPdn(:) = R*T/V - ArVn(:)
   end subroutine pressure

   subroutine fugacity_pt(eos, &
      n, P, T, V, root_type, lnPhi, dlnPhidP, dlnPhidT, dlnPhidn, dPdV, dPdT, dPdn &
      )
      !! Calculate logarithm of fugacity, given pressure and temperature.
      !!
      !! This routine will obtain the desired volume root at the specified
      !! pressure and calculate fugacity at that point.
      use iso_fortran_env, only: error_unit
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Mixture mole numbers
      character(len=*), intent(in) :: root_type !! Type of root desired ["liquid", "vapor", "stable"]
      real(pr), intent(in) :: P    !! Pressure [bar]
      real(pr), intent(in) :: T    !! Temperature [K]

      real(pr), optional, intent(out) :: lnPhi(size(n)) !! \(\ln(phi)\) vector
      real(pr), optional, intent(out) :: V !! Volume [L]
      real(pr), optional, intent(out) :: dlnPhidT(size(n)) !! ln(phi) Temp derivative
      real(pr), optional, intent(out) :: dlnPhidP(size(n)) !! ln(phi) Presssure derivative
      real(pr), optional, intent(out) :: dlnPhidn(size(n), size(n)) !! ln(phi) compositional derivative
      real(pr), optional, intent(out) :: dPdV !! \(\frac{dP}{dV}\)
      real(pr), optional, intent(out) :: dPdT !! \(\frac{dP}{dT}\)
      real(pr), optional, intent(out) :: dPdn(size(n)) !! \(\frac{dP}{dn_i}\)

      real(pr) :: V_in, P_in

      call eos%volume(n, P=P, T=T, V=V_in, root_type=root_type)
      call eos%lnphi_vt(&
         n, V=V_in, T=T, &
         P=P_in, lnPhi=lnPhi, &
         dlnPhidP=dlnPhidP, dlnPhidT=dlnPhidT, dlnPhidn=dlnPhidn, &
         dPdV=dPdV, dPdT=dPdT, dPdn=dPdn &
         )

      if(present(V)) V = V_in

      ! Check if the calculated pressure is the same as the input pressure.
      if(abs(P_in - P) > 1e-2) then
         write(error_unit, *) "WARN: possible bad root solving: ", P_in, P
      end if
   end subroutine fugacity_pt

   subroutine fugacity_vt(eos, &
      n, V, T, P, lnPhi, &
      dlnPhidP, dlnPhidT, dlnPhidn, &
      dPdV, dPdT, dPdn &
      )
      !! Calculate fugacity coefficent given volume and temperature.
      class(ArModel) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Mixture mole numbers
      real(pr), intent(in) :: V !! Volume [L]
      real(pr), intent(in) :: T !! Temperature [K]

      real(pr), optional, intent(out) :: P !! Pressure [bar]
      real(pr), optional, intent(out) :: lnPhi(size(n)) !! \(\ln(\phi_i)\) vector
      real(pr), optional, intent(out) :: dlnPhidT(size(n)) !! \(ln(phi_i)\) Temp derivative
      real(pr), optional, intent(out) :: dlnPhidP(size(n)) !! \(ln(phi_i)\) Presssure derivative
      real(pr), optional, intent(out) :: dlnPhidn(size(n), size(n)) !! \(ln(phi_i)\) compositional derivative
      real(pr), optional, intent(out) :: dPdV !! \(\frac{dP}{dV}\)
      real(pr), optional, intent(out) :: dPdT !! \(\frac{dP}{dT}\)
      real(pr), optional, intent(out) :: dPdn(:) !! \(\frac{dP}{dn_i}\)

      real(pr) :: Ar, ArTV, ArV, ArV2
      real(pr), dimension(size(n)) :: Arn, ArVn, ArTn
      real(pr) :: Arn2(size(n), size(n))

      real(pr) :: dPdV_in, dPdT_in, dPdn_in(size(n))
      real(pr) :: P_in

      real(pr) :: RT, Z

      real(pr) :: totn
      integer :: nc, i, j

      totn = sum(n)
      nc = size(n)

      RT = R*T

      if (present(lnPhi) .and. .not. (&
         present(dlnPhidn) &
         .or. present(dlnPhidP) &
         .or. present(dlnPhidT) &
         )) then
         call eos%residual_helmholtz(n, v, t, Arn=Arn, ArV=ArV)
         P_in = totn*RT/V - ArV
         Z = P_in*V/(totn*RT)
         lnPhi(:) = Arn(:)/RT - log(Z)
         if (present(P)) P = P_in
         return
      else if (present(dlnPhidn)) then
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

      P_in = totn*RT/V - ArV

      Z = P_in*V/(totn*RT)
      if (present(P)) P = P_in

      dPdV_in = -ArV2 - RT*totn/V**2
      dPdT_in = -ArTV + totn*R/V
      dPdn_in = RT/V - ArVn

      if (present(lnPhi)) lnPhi = Arn(:)/RT - log(Z)
      if (present(dlnPhidP)) then
         dlnPhidP(:) = -dPdn_in(:)/dPdV_in/RT - 1._pr/P_in
      end if
      if (present(dlnPhidT)) then
         dlnPhidT(:) = (ArTn(:) - Arn(:)/T)/RT + dPdn_in(:)*dPdT_in/dPdV_in/RT + 1._pr/T
      end if

      if (present(dlnPhidn)) then
         do i = 1, nc
            do j = i, nc
               dlnPhidn(i, j) = 1._pr/totn + (Arn2(i, j) + dPdn_in(i)*dPdn_in(j)/dPdV_in)/RT
               dlnPhidn(j, i) = dlnPhidn(i, j)
            end do
         end do
      end if

      if (present(dPdV)) dPdV = dPdV_in
      if (present(dPdT)) dPdT = dPdT_in
      if (present(dPdn)) dPdn = dPdn_in
   end subroutine fugacity_vt

   subroutine enthalpy_residual_vt(eos, n, V, T, Hr, HrT, HrV, Hrn)
      !! Calculate residual enthalpy given volume and temperature.
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: t !! Temperature [K]
      real(pr), intent(in) :: v !! Volume [L]
      real(pr), intent(out) :: Hr !! Residual enthalpy [bar L]
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

   subroutine gibbs_residual_VT(eos, n, V, T, Gr, GrT, GrV, Grn)
      !! Calculate residual Gibbs energy given volume and temperature.
      use yaeos__constants, only: R
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: V !! Volume [L]
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(out) :: Gr !! Gibbs energy [bar L]
      real(pr), optional, intent(out) :: GrT !! \(\frac{dG^r}}{dT}\)
      real(pr), optional, intent(out) :: GrV !! \(\frac{dG^r}}{dV}\)
      real(pr), optional, intent(out) :: Grn(size(n)) !! \(\frac{dG^r}}{dn}\)

      real(pr) :: Ar, ArV, ArT, Arn(size(n))
      real(pr) :: p, dPdV, dPdT, dPdn(size(n)), z, totn

      totn = sum(n)
      call pressure(eos, n, V, T, P, dPdV=dPdV, dPdT=dPdT, dPdn=dPdn)
      z = P*V/(totn*R*T)

      call eos%residual_helmholtz(n, v, t, Ar=Ar, ArV=ArV, ArT=ArT, Arn=Arn)

      Gr = Ar + P*V - totn*R*T

      if (present(GrT)) GrT = ArT + V*dPdT - totn*R
      if (present(GrV)) GrV = ArV + V*dPdV + P
      if (present(GrN)) GrN(:) = Arn(:) + V*dPdn(:) - R*T
   end subroutine gibbs_residual_VT

   subroutine entropy_residual_vt(eos, n, V, T, Sr, SrT, SrV, Srn)
      !! Calculate residual entropy given volume and temperature.
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: V !! Volume [L]
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(out) :: Sr !! Entropy [bar L / K]
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

   subroutine Cv_residual_vt(eos, n, V, T, Cv)
      !! Calculate residual heat capacity volume constant given v and t.
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: V !! Volume [L]
      real(pr), intent(out) :: Cv !! heat capacity v constant [bar L / K]

      real(pr) :: Ar, ArT2

      call eos%residual_helmholtz(n, V, T, Ar=Ar, ArT2=ArT2)

      Cv = -T*ArT2
   end subroutine Cv_residual_vt

   subroutine Cp_residual_vt(eos, n, V, T, Cp)
      !! Calculate residual heat capacity pressure constant given v and t.
      use yaeos__constants, only: R
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: V !! Volume [L]
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(out) :: Cp !! heat capacity p constant [bar L / K]

      real(pr) :: Ar, ArT2, Cv, p, dPdT, dPdV, totn

      totn = sum(n)

      call eos%residual_helmholtz(n, V, T, Ar=Ar, ArT2=ArT2)

      call Cv_residual_vt(eos, n, V, T, Cv)

      call pressure(eos, n, V, T, P, dPdV=dPdV, dPdT=dPdT)

      Cp = Cv - T*dPdT**2/dPdV - totn*R
   end subroutine Cp_residual_vt
end module yaeos__models_ar
