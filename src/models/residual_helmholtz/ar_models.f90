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

   private

   public :: ArModel, size, volume

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
      procedure :: lnphi_vt
      procedure :: lnphi_pt
      procedure :: lnfug_vt
      procedure :: pressure
      procedure :: volume
      procedure :: enthalpy_residual_vt
      procedure :: gibbs_residual_vt
      procedure :: entropy_residual_vt
      procedure :: internal_energy_residual_vt
      procedure :: Cv_residual_vt
      procedure :: Cp_residual_vt
      procedure :: helmholtz_residual_pt
      procedure :: enthalpy_residual_pt
      procedure :: gibbs_residual_pt
      procedure :: entropy_residual_pt
      procedure :: internal_energy_residual_pt
      procedure :: Cv_residual_pt
      procedure :: Cp_residual_pt
      procedure :: ln_activity_coefficient
      procedure :: gibbs_excess
      procedure :: enthalpy_excess
      procedure :: volume_excess
      procedure :: entropy_excess
      procedure :: Psat_pure
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
      !! Volume solver routine for residual Helmholtz models.
      !!
      !! Solves volume roots using newton method. Given pressure and
      !! temperature.
      !!
      !! # Description
      !! This subroutine solves the volume using a newton method. The variable
      !! `root_type` is used to specify the desired root to solve. The options
      !! are: `["liquid", "vapor", "stable"]`
      !!
      !! # Examples
      !!
      !! ```fortran
      !! eos = PengRobinson76(Tc, Pc, w)
      !!
      !! n = [1.0_pr, 1.0_pr]
      !! T = 300.0_pr
      !! P = 1.0_pr
      !!
      !! call eos%volume(n, P, T, V, root_type="liquid")
      !! call eos%volume(n, P, T, V, root_type="vapor")
      !! call eos%volume(n, P, T, V, root_type="stable")
      !! ```
      use yaeos__math, only: newton

      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(out) :: V !! Volume [L]
      character(len=*), intent(in) :: root_type
      !! Desired root-type to solve. Options are:
      !! `["liquid", "vapor", "stable"]`

      integer :: max_iters=30
      real(pr) :: tol=1e-8

      real(pr) :: totnRT, GrL, GrV, Gr
      real(pr) :: Vliq, Vvap
      logical :: failed

      GrL = HUGE(GrL)
      GrV = HUGE(GrV)

      totnRT = sum(n) * R * T
      select case(root_type)
       case("liquid")
         Vliq = eos%get_v0(n, P, T) * 1.001_pr
         call newton(foo, Vliq, tol=tol, max_iters=max_iters, failed=failed)
         GrL = Gr
       case("vapor")
         ! Vvap = R * T / P
         Vvap = R * T / P
         call newton(foo, Vvap, tol=tol, max_iters=max_iters, failed=failed)
         GrV = Gr
       case("stable")
         Vliq = eos%get_v0(n, P, T)*1.00001_pr
         call newton(foo, Vliq, tol=tol, max_iters=max_iters, failed=failed)
         GrL = Gr

         Vvap = R * T / P
         call newton(foo, Vvap, tol=tol, max_iters=max_iters, failed=failed)
         GrV = Gr
      end select

      if (GrL < GrV) then
         V = Vliq
      else
         V = Vvap
      end if

      if (failed) V = -1

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
      !! Calculate pressure.
      !!
      !! Calculate pressure using residual helmholtz models.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! eos = PengRobinson76(Tc, Pc, w)
      !!
      !! n = [1.0_pr, 1.0_pr]
      !! T = 300.0_pr
      !! V = 1.0_pr
      !!
      !! call eos%pressure(n, V, T, P, dPdV=dPdV, dPdT=dPdT, dPdn=dPdn)
      !! ```
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: V !! Volume [L]
      real(pr), intent(out) :: P !! Pressure [bar]
      real(pr), optional, intent(out) :: dPdV !! \(\frac{dP}}{dV}\)
      real(pr), optional, intent(out) :: dPdT !! \(\frac{dP}}{dT}\)
      real(pr), optional, intent(out) :: dPdn(:) !! \(\frac{dP}}{dn_i}\)

      real(pr) :: totn

      real(pr) :: Ar, ArV, ArV2, ArTV, ArVn(size(eos))

      totn = sum(n)

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

   subroutine lnphi_pt(eos, &
      n, P, T, V, root_type, lnPhi, dlnPhidP, dlnPhidT, dlnPhidn, &
      dPdV, dPdT, dPdn, lnPhiP &
      )
      !! Calculate natural logarithm of fugacity given pressure and temperature.
      !!
      !! Calculate the natural logarithm of the fugacity coefficient and its
      !! derivatives given pressure and temperature. This routine will obtain
      !! the desired volume root at the specified pressure and calculate
      !! fugacity at that point.The routine gives the possibility to calculate
      !! the pressure derivatives and volume.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! eos = PengRobinson76(Tc, Pc, w)
      !!
      !! n = [1.0_pr, 1.0_pr]
      !! T = 300.0_pr
      !! V = 1.0_pr
      !!
      !! call eos%lnphi_pt(&
      !!    n, V, T, lnPhi=lnPhi, &
      !!    dlnPhidP=dlnPhidP, dlnPhidT=dlnPhidT, dlnPhidn=dlnPhidn &
      !!    )
      !! ```
      use iso_fortran_env, only: error_unit

      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Mixture mole numbers
      character(len=*), intent(in) :: root_type
      !! Type of root desired ["liquid", "vapor", "stable"]
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]

      real(pr), optional, intent(out) :: lnPhi(size(n)) !! \(\ln(phi)\) vector
      real(pr), optional, intent(out) :: V !! Volume [L]
      real(pr), optional, intent(out) :: dlnPhidT(size(n)) !! ln(phi) Temperature derivative
      real(pr), optional, intent(out) :: dlnPhidP(size(n)) !! ln(phi) Presssure derivative
      real(pr), optional, intent(out) :: dlnPhidn(size(n), size(n)) !! ln(phi) compositional derivative
      real(pr), optional, intent(out) :: dPdV !! \(\frac{dP}{dV}\)
      real(pr), optional, intent(out) :: dPdT !! \(\frac{dP}{dT}\)
      real(pr), optional, intent(out) :: dPdn(size(n)) !! \(\frac{dP}{dn_i}\)
      real(pr), optional, intent(out) :: lnPhiP(:)
      !! \(\ln(\phi_i)P \). It is useful to calculate fugacity coefficients
      !! at negatives pressures.

      real(pr) :: V_in, P_in

      call eos%volume(n, P=P, T=T, V=V_in, root_type=root_type)

      call eos%lnphi_vt(&
         n, V=V_in, T=T, &
         P=P_in, lnPhi=lnPhi, &
         dlnPhidP=dlnPhidP, dlnPhidT=dlnPhidT, dlnPhidn=dlnPhidn, &
         dPdV=dPdV, dPdT=dPdT, dPdn=dPdn, lnPhiP=lnPhiP &
         )

      if(present(V)) V = V_in

      ! Check if the calculated pressure is the same as the input pressure.
      ! if(abs(P_in - P) > 1e-2) then
      !    write(error_unit, *) "WARN: possible bad root solving: ", P_in, P
      ! end if
   end subroutine lnphi_pt

   subroutine lnphi_vt(eos, &
      n, V, T, P, lnPhi, &
      dlnPhidP, dlnPhidT, dlnPhidn, &
      dPdV, dPdT, dPdn, lnPhiP &
      )
      !! Calculate natural logarithm of fugacity coefficent.
      !!
      !! Calculate the natural logarithm of the fugacity coefficient and its
      !! derivatives given volume and temperature. The routine gives the
      !! possibility to calculate the pressure and it's derivatives.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! eos = PengRobinson76(Tc, Pc, w)
      !!
      !! n = [1.0_pr, 1.0_pr]
      !! T = 300.0_pr
      !! V = 1.0_pr
      !!
      !! call eos%lnphi_vt(&
      !!    n, V, T, lnPhi=lnPhi, &
      !!    dlnPhidP=dlnPhidP, dlnPhidT=dlnPhidT, dlnPhidn=dlnPhidn &
      !!    )
      !! ```
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
      real(pr), optional, intent(out) :: lnPhiP(:)
      !! \(\ln(\phi_i)P \). It is useful to calculate fugacity coefficients
      !! at negatives pressures.

      real(pr) :: Ar, ArTV, ArV, ArV2
      real(pr), dimension(size(n)) :: Arn, ArVn, ArTn
      real(pr) :: Arn2(size(n), size(n))

      real(pr) :: dPdV_in, dPdT_in, dPdn_in(size(n))
      real(pr) :: P_in

      real(pr) :: RT, Z, Z_P

      real(pr) :: totn
      integer :: nc, i, j

      totn = sum(n)
      nc = size(n)

      RT = R*T

      if ((present(lnPhi) .or. present(lnPhiP)) .and. .not. (&
         present(dlnPhidn) &
         .or. present(dlnPhidP) &
         .or. present(dlnPhidT) &
         )) then
         call eos%residual_helmholtz(n, v, t, Arn=Arn, ArV=ArV)

         P_in = totn*RT/V - ArV
         if (present(P)) P = P_in

         if (present(lnPhi)) then
            Z = P_in*V/(totn*RT)
            lnPhi(:) = Arn(:)/RT - log(Z)
         end if

         if (present(lnPhiP)) then
            Z = V/(totn*RT)
            lnPhiP(:) = Arn(:)/RT - log(Z)
         end if
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
      if (present(lnPhiP)) then
         ! Avoiding P in Z makes it possible to calculate fugacities at
         ! negative pressures
         Z_P = V/(totn*RT)
         lnPhiP = Arn(:)/RT - log(Z_P)
      end if

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
   end subroutine lnphi_vt

   subroutine lnfug_vt(eos, &
      n, V, T, P, lnf, &
      dlnfdV, dlnfdT, dlnfdn, &
      dPdV, dPdT, dPdn &
      )
      !! Calculate natural logarithm of fugacity given volume and temperature.
      !!
      !! Calculate the natural logarithm of the fugacity and its derivatives.
      !! The routine gives the possibility to calculate the pressure and it's
      !! derivatives.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! eos = PengRobinson76(Tc, Pc, w)
      !!
      !! n = [1.0_pr, 1.0_pr]
      !! T = 300.0_pr
      !! V = 1.0_pr
      !!
      !! call eos%lnfug_vt(&
      !!    n, V, T, lnf=lnf, &
      !!    dlnfdV=dlnfdV, dlnfdT=dlnfdT, dlnfdn=dlnfdn &
      !!    )
      !! ```
      class(ArModel) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Mixture mole numbers
      real(pr), intent(in) :: V !! Volume [L]
      real(pr), intent(in) :: T !! Temperature [K]

      real(pr), optional, intent(out) :: P !! Pressure [bar]
      real(pr), optional, intent(out) :: lnf(size(n)) !! \(\ln(\f_i)\) vector
      real(pr), optional, intent(out) :: dlnfdT(size(n)) !! \(ln(f_i)\) Temp derivative
      real(pr), optional, intent(out) :: dlnfdV(size(n)) !! \(ln(f_i)\) Volume derivative
      real(pr), optional, intent(out) :: dlnfdn(size(n), size(n)) !! \(ln(f_i)\) compositional derivative
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

      if (present(lnf) .and. .not. (&
         present(dlnfdn) &
         .or. present(dlnfdV) &
         .or. present(dlnfdT) &
         )) then
         call eos%residual_helmholtz(n, v, t, Arn=Arn, ArV=ArV)

         P_in = totn*RT/V - ArV

         lnf = 0.0_pr

         where (n /= 0)
            lnf = log(n/totn) + Arn/RT - log(V/(totn*RT))
         endwhere

         if (present(P)) P = P_in

         return
      else if (present(dlnfdn)) then
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

      if (present(lnf)) then
         lnf = 0.0_pr

         where (n /= 0)
            lnf = log(n/totn) + Arn/RT - log(V/(totn*RT))
         endwhere
      end if

      if (present(dlnfdV)) then
         dlnfdV = -dPdn_in/RT
      end if

      if (present(dlnfdT)) then
         dlnfdT = (ArTn - Arn/T)/RT + 1._pr/T
      end if

      if (present(dlnfdn)) then
         do i = 1, nc
            do j=1,nc
               dlnfdn(i, j) = Arn2(i, j)/RT
            end do
            dlnfdn(i, i) = dlnfdn(i, i) + 1/n(i)
         end do
      end if

      if (present(dPdV)) dPdV = dPdV_in
      if (present(dPdT)) dPdT = dPdT_in
      if (present(dPdn)) dPdn = dPdn_in
   end subroutine lnfug_vt

   ! ==========================================================================
   ! VT thermoprops
   ! --------------------------------------------------------------------------
   subroutine enthalpy_residual_vt(eos, n, V, T, Hr, HrV, HrT, Hrn)
      !! Calculate residual enthalpy given volume and temperature.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! eos = PengRobinson76(Tc, Pc, w)
      !!
      !! n = [1.0_pr, 1.0_pr]
      !! T = 300.0_pr
      !! V = 1.0_pr
      !!
      !! call eos%enthalpy_residual_vt(&
      !!    n, V, T, Hr=Hr, HrV=HrV, HrT=HrT, Hrn=Hrn &
      !!    )
      !! ```
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: V !! Volume [L]
      real(pr), optional, intent(out) :: Hr !! Residual enthalpy [bar L]
      real(pr), optional, intent(out) :: HrV !! \(\frac{dH^r}}{dV}\)
      real(pr), optional, intent(out) :: HrT !! \(\frac{dH^r}}{dT}\)
      real(pr), optional, intent(out) :: Hrn(size(n)) !! \(\frac{dH^r}}{dn}\)

      real(pr) :: Ar, ArV, ArT, Arn(size(n))
      real(pr) :: ArV2, ArT2, ArTV, ArVn(size(n)), ArTn(size(n))

      call eos%residual_helmholtz(&
         n, V, T, Ar=Ar, ArV=ArV, ArT=ArT, ArTV=ArTV, &
         ArV2=ArV2, ArT2=ArT2, Arn=Arn, ArVn=ArVn, ArTn=ArTn &
         )

      if (present(Hr)) Hr = Ar - T*ArT - V*ArV
      if (present(HrT)) HrT = -T*ArT2 - V*ArTV
      if (present(HrV)) HrV = -T*ArTV - V*ArV2
      if (present(Hrn)) Hrn(:) = Arn(:) - T*ArTn(:) - V*ArVn(:)
   end subroutine enthalpy_residual_vt

   subroutine gibbs_residual_vt(eos, n, V, T, Gr, GrV, GrT, Grn)
      !! Calculate residual Gibbs energy given volume and temperature.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! eos = PengRobinson76(Tc, Pc, w)
      !!
      !! n = [1.0_pr, 1.0_pr]
      !! T = 300.0_pr
      !! V = 1.0_pr
      !!
      !! call eos%gibbs_residual_vt(&
      !!    n, V, T, Gr=Gr, GrV=GrV, GrT=GrT, Grn=Grn &
      !!    )
      !! ```
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: V !! Volume [L]
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), optional, intent(out) :: Gr !! Gibbs energy [bar L]
      real(pr), optional, intent(out) :: GrV !! \(\frac{dG^r}}{dV}\)
      real(pr), optional, intent(out) :: GrT !! \(\frac{dG^r}}{dT}\)
      real(pr), optional, intent(out) :: Grn(size(n)) !! \(\frac{dG^r}}{dn}\)

      real(pr) :: Ar, ArV, ArT, Arn(size(n)), ArTV, ArV2, ArVn(size(n))

      call eos%residual_helmholtz(&
         n, v, t, Ar=Ar, ArV=ArV, &
         ArT=ArT, Arn=Arn, ArTV=ArTV, ArV2=ArV2, ArVn=ArVn &
         )

      if (present(Gr)) Gr = Ar - V*ArV
      if (present(GrT)) GrT = ArT - V*ArTV
      if (present(GrV)) GrV = -V*ArV2
      if (present(Grn)) Grn(:) = Arn(:) - V*ArVn(:)
   end subroutine gibbs_residual_vt

   subroutine entropy_residual_vt(eos, n, V, T, Sr, SrV, SrT, Srn)
      !! Calculate residual entropy given volume and temperature.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! eos = PengRobinson76(Tc, Pc, w)
      !!
      !! n = [1.0_pr, 1.0_pr]
      !! T = 300.0_pr
      !! V = 1.0_pr
      !!
      !! call eos%entropy_residual_vt(&
      !!    n, V, T, Sr=Sr, SrV=SrV, SrT=SrT, Srn=Srn &
      !!    )
      !! ```
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: V !! Volume [L]
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), optional, intent(out) :: Sr !! Entropy [bar L / K]
      real(pr), optional, intent(out) :: SrV !! \(\frac{dS^r}}{dV}\)
      real(pr), optional, intent(out) :: SrT !! \(\frac{dS^r}}{dT}\)
      real(pr), optional, intent(out) :: Srn(size(n)) !! \(\frac{dS^r}}{dn}\)

      real(pr) :: Ar, ArT, ArT2, ArTV, ArTn(size(n))

      call eos%residual_helmholtz(&
         n, v, t, Ar=Ar, ArT=ArT, ArTV=ArTV, ArT2=ArT2, ArTn=ArTn &
         )

      if (present(Sr)) Sr = -ArT
      if (present(SrT)) SrT = -ArT2
      if (present(SrV)) SrV = -ArTV
      if (present(SrN)) Srn = -ArTn
   end subroutine entropy_residual_vt

   subroutine internal_energy_residual_vt(eos, n, V, T, Ur, UrV, UrT, Urn)
      !! Calculate residual internal energy given volume and temperature.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! eos = PengRobinson76(Tc, Pc, w)
      !!
      !! n = [1.0_pr, 1.0_pr]
      !! T = 300.0_pr
      !! V = 1.0_pr
      !!
      !! call eos%internal_energy_residual_vt(&
      !!    n, V, T, Ur=Ur, UrV=UrV, UrT=UrT, Urn=Urn &
      !!    )
      !! ```
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: V !! Volume [L]
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), optional, intent(out) :: Ur !! Internal energy [bar L]
      real(pr), optional, intent(out) :: UrV !! \(\frac{dU^r}}{dV}\)
      real(pr), optional, intent(out) :: UrT !! \(\frac{dU^r}}{dT}\)
      real(pr), optional, intent(out) :: Urn(size(n)) !! \(\frac{dU^r}}{dn}\)

      real(pr) :: Ar, ArV, ArT, Arn(size(n)), ArT2, ArTV, ArTn(size(n))

      call eos%residual_helmholtz(&
         n, v, t, Ar=Ar, ArV=ArV, ArT=ArT, Arn=Arn, &
         ArTV=ArTV, ArT2=ArT2, ArTn=ArTn &
         )

      if (present(Ur)) Ur = Ar - T*ArT
      if (present(UrT)) UrT = -T*ArT2
      if (present(UrV)) UrV = ArV - T*ArTV
      if (present(Urn)) Urn = Arn - T*ArTn
   end subroutine internal_energy_residual_vt

   subroutine Cv_residual_vt(eos, n, V, T, Cv)
      !! Calculate residual heat capacity volume constant given V and T.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! eos = PengRobinson76(Tc, Pc, w)
      !!
      !! n = [1.0_pr, 1.0_pr]
      !! T = 300.0_pr
      !! V = 1.0_pr
      !!
      !! call eos%Cv_residual_vt(n, V, T, Cv=Cv)
      !! ```
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: V !! Volume [L]
      real(pr), intent(out) :: Cv !! heat capacity V constant [bar L / K]

      real(pr) :: ArT2

      call eos%residual_helmholtz(n, V, T, ArT2=ArT2)

      Cv = -T*ArT2
   end subroutine Cv_residual_vt

   subroutine Cp_residual_vt(eos, n, V, T, Cp)
      !! Calculate residual heat capacity pressure constant given V and T.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! eos = PengRobinson76(Tc, Pc, w)
      !!
      !! n = [1.0_pr, 1.0_pr]
      !! T = 300.0_pr
      !! V = 1.0_pr
      !!
      !! call eos%Cp_residual_vt(n, V, T, Cp=Cp)
      !! ```
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: V !! Volume [L]
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(out) :: Cp !! heat capacity P constant [bar L / K]

      real(pr) :: nt
      real(pr) :: ArT2, ArTV, ArV2

      nt = sum(n)

      call eos%residual_helmholtz(n, V, T, ArTV=ArTV, ArV2=ArV2, ArT2=ArT2)

      Cp = &
         -T * ArT2 &
         -T * (-ArTV + nt *R / V)**2 / (-ArV2 - nt*R*T / V**2) &
         - nt * R
   end subroutine Cp_residual_vt

   ! ==========================================================================
   ! PT thermoprops
   ! --------------------------------------------------------------------------
   subroutine helmholtz_residual_pt(eos, n, P, T, root_type, Ar, ArP, ArT, Arn)
      !! Calculate residual Helmholtz energy given pressure and temperature.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! eos = PengRobinson76(Tc, Pc, w)
      !!
      !! n = [1.0_pr, 1.0_pr]
      !! T = 300.0_pr
      !! P = 1.0_pr
      !!
      !! call eos%helmholtz_residual_pt(&
      !!    n, P, T, root_type="stable", Ar=Ar, ArP=ArP, ArT=ArT, Arn=Arn &
      !!    )
      !! ```
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]
      character(len=*), intent(in) :: root_type
      !! Desired root-type to solve. Options are:
      !! `["liquid", "vapor", "stable"]`
      real(pr), optional, intent(out) :: Ar
      !! Residual Helmholtz energy [bar L]
      real(pr), optional, intent(out) :: ArP !! \(\frac{dA^r}}{dP}\)
      real(pr), optional, intent(out) :: ArT !! \(\frac{dA^r}}{dT}\)
      real(pr), optional, intent(out) :: Arn(size(n)) !! \(\frac{dA^r}}{dn}\)

      ! Helmholtz at VT
      real(pr) :: Ar_v, ArV_v, ArT_v, Arn_v(size(n))
      real(pr) :: dPdV, dPdT, dPdn(size(n))
      real(pr) :: dVdP, dVdT, dVdn(size(n))
      real(pr) :: Z, nt, V, P_dummy

      logical :: dp, dt, dn, present_derivs

      dp = present(ArP)
      dt = present(ArT)
      dn = present(Arn)

      present_derivs = dp .or. dt .or. dn

      if (present_derivs) then
         call eos%volume(n=n, P=P, T=T, V=V, root_type=root_type)
         call eos%pressure(&
            n=n, V=V, T=T, P=P_dummy, dPdV=dPdV, dPdT=dPdT, dPdn=dPdn &
            )
         call eos%residual_helmholtz(&
            n=n, V=V, T=T, Ar=Ar_v, ArV=ArV_v, ArT=ArT_v, Arn=Arn_v &
            )

         dVdP = 1 / dPdV
         dVdT = -dPdT / dPdV
         dVdn = -dPdn / dPdV
      else
         call eos%volume(n=n, P=P, T=T, V=V, root_type=root_type)
         call eos%residual_helmholtz(n=n, V=V, T=T, Ar=Ar_v)
      end if

      nt = sum(n)
      Z = P * V / (nt * R * T)

      if (present(Ar)) Ar = Ar_v - nt * R * T * log(Z)

      if (present(ArP)) ArP = ArV_v / dPdV - nt*R*T*(1.0_pr/P + dVdP/V)

      if (present(ArT)) ArT = &
         ArT_v + ArV_v * dVdT - nt*R*log(Z) - nt*R*T*(dVdT / V - 1.0_pr/T)

      if (present(Arn)) Arn = &
         Arn_v + ArV_v * dVdn - R*T*log(Z) - nt*R*T*(dVdn / V - 1.0_pr/nt)
   end subroutine helmholtz_residual_pt

   subroutine enthalpy_residual_pt(eos, n, P, T, root_type, Hr, HrP, HrT, Hrn)
      !! Calculate residual enthalpy given pressure and temperature.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! eos = PengRobinson76(Tc, Pc, w)
      !!
      !! n = [1.0_pr, 1.0_pr]
      !! T = 300.0_pr
      !! P = 1.0_pr
      !!
      !! call eos%enthalpy_residual_pt(&
      !!    n, P, T, root_type="stable", Hr=Hr, HrP=HrP, HrT=HrT, Hrn=Hrn &
      !!    )
      !! ```
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]
      character(len=*), intent(in) :: root_type
      !! Desired root-type to solve. Options are:
      !! `["liquid", "vapor", "stable"]`
      real(pr), optional, intent(out) :: Hr !! Residual enthalpy [bar L]
      real(pr), optional, intent(out) :: HrP !! \(\frac{dH^r}}{dP}\)
      real(pr), optional, intent(out) :: HrT !! \(\frac{dH^r}}{dT}\)
      real(pr), optional, intent(out) :: Hrn(size(n)) !! \(\frac{dH^r}}{dn}\)

      real(pr) :: Hr_v, HrT_v, HrV_v, Hrn_v(size(n)) ! residual entralpy vt
      real(pr) :: V, dPdV, dPdT, dPdn(size(n)), P_dummy
      real(pr) :: dVdP, dVdT, dVdn(size(n))

      logical :: dp, dT, dn, derivs_present

      dp = present(HrP)
      dt = present(HrT)
      dn = present(Hrn)

      derivs_present = dp .or. dt .or. dn

      if (derivs_present) then
         call eos%volume(n=n, P=P, T=T, V=V, root_type=root_type)
         call eos%pressure(&
            n=n, V=V, T=T, P=P_dummy, dPdV=dPdV, dPdT=dPdT, dPdn=dPdn &
            )
         call eos%enthalpy_residual_vt(&
            n=n, V=V, T=T, Hr=Hr_v, HrV=HrV_v, HrT=HrT_v, Hrn=Hrn_v &
            )

         dVdP = 1 / dPdV
         dVdT = -dPdT / dPdV
         dVdn = -dPdn / dPdV
      else
         call eos%volume(n=n, P=P, T=T, V=V, root_type=root_type)
         call eos%enthalpy_residual_vt(n=n, V=V, T=T, Hr=Hr_v)
      end if

      if (present(Hr)) Hr = Hr_v
      if (present(HrP)) HrP = HrV_v * dVdP
      if (present(HrT)) HrT = HrT_v + HrV_v * dVdT
      if (present(Hrn)) Hrn = Hrn_v + HrV_v * dVdn
   end subroutine enthalpy_residual_pt

   subroutine gibbs_residual_pt(eos, n, P, T, root_type, Gr, GrP, GrT, Grn)
      !! Calculate residual Gibbs energy given pressure and temperature.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! eos = PengRobinson76(Tc, Pc, w)
      !!
      !! n = [1.0_pr, 1.0_pr]
      !! T = 300.0_pr
      !! P = 1.0_pr
      !!
      !! call eos%gibbs_residual_pt(&
      !!    n, P, T, root_type="stable", Gr=Gr, GrP=GrP, GrT=GrT, Grn=Grn &
      !!    )
      !! ```
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]
      character(len=*), intent(in) :: root_type
      !! Desired root-type to solve. Options are:
      !! `["liquid", "vapor", "stable"]`
      real(pr), optional, intent(out) :: Gr !! Residual Gibbs energy [bar L]
      real(pr), optional, intent(out) :: GrP !! \(\frac{dG^r}}{dP}\)
      real(pr), optional, intent(out) :: GrT !! \(\frac{dG^r}}{dT}\)
      real(pr), optional, intent(out) :: Grn(size(n)) !! \(\frac{dG^r}}{dn}\)

      real(pr) :: Gr_v, GrT_v, GrV_v, Grn_v(size(n)) ! residual Gibbs energy vt
      real(pr) :: Z, nt, V, dPdV, dPdT, dPdn(size(n)), P_dummy
      real(pr) :: dVdP, dVdT, dVdn(size(n))

      logical :: dp, dt, dn, present_derivs

      dp = present(GrP)
      dt = present(GrT)
      dn = present(Grn)

      present_derivs = dp .or. dt .or. dn

      if (present_derivs) then
         call eos%volume(n=n, P=P, T=T, V=V, root_type=root_type)
         call eos%pressure(&
            n=n, V=V, T=T, P=P_dummy, dPdV=dPdV, dPdT=dPdT, dPdn=dPdn &
            )
         call eos%gibbs_residual_vt(&
            n=n, V=V, T=T, Gr=Gr_v, GrV=GrV_v, GrT=GrT_v, Grn=Grn_v &
            )

         dVdP = 1 / dPdV
         dVdT = -dPdT / dPdV
         dVdn = -dPdn / dPdV
      else
         call eos%volume(n=n, P=P, T=T, V=V, root_type=root_type)
         call eos%gibbs_residual_vt(n=n, V=V, T=T, Gr=Gr_v)
      end if

      nt = sum(n)
      Z = P * V / (nt * R * T)

      if (present(Gr)) Gr = Gr_v - nt * R * T * log(Z)

      if (present(GrP)) GrP = GrV_v / dPdV - nt*R*T*(1.0_pr/P + dVdP/V)

      if (present(GrT)) GrT = &
         GrT_v + GrV_v * dVdT - nt*R*log(Z) - nt*R*T*(dVdT / V - 1.0_pr/T)

      if (present(Grn)) Grn = &
         Grn_v + GrV_v * dVdn - R*T*log(Z) - nt*R*T*(dVdn / V - 1.0_pr/nt)
   end subroutine gibbs_residual_pt

   subroutine entropy_residual_pt(eos, n, P, T, root_type, Sr, SrP, SrT, Srn)
      !! Calculate residual entropy given pressure and temperature.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! eos = PengRobinson76(Tc, Pc, w)
      !!
      !! n = [1.0_pr, 1.0_pr]
      !! T = 300.0_pr
      !! P = 1.0_pr
      !!
      !! call eos%entropy_residual_pt(&
      !!    n, P, T, root_type="stable", Sr=Sr, SrP=SrP, SrT=SrT, Srn=Srn &
      !!    )
      !! ```
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]
      character(len=*), intent(in) :: root_type
      !! Desired root-type to solve. Options are:
      !! `["liquid", "vapor", "stable"]`
      real(pr), optional, intent(out) :: Sr !! Residual entropy [bar L K^-1]
      real(pr), optional, intent(out) :: SrP !! \(\frac{dS^r}}{dP}\)
      real(pr), optional, intent(out) :: SrT !! \(\frac{dS^r}}{dT}\)
      real(pr), optional, intent(out) :: Srn(size(n)) !! \(\frac{dS^r}}{dn}\)

      real(pr) :: Sr_v, SrT_v, SrV_v, Srn_v(size(n)) ! residual entropy vt
      real(pr) :: dPdV, dPdT, dPdn(size(n))
      real(pr) :: dVdT, dVdP, dVdn(size(n))
      real(pr) :: Z, nt, V, P_dummy

      logical :: dp, dt, dn, derivs_present

      dp = present(SrP)
      dt = present(SrT)
      dn = present(Srn)
      derivs_present = dp .or. dt .or. dn

      if (derivs_present) then
         call eos%volume(n=n, P=P, T=T, V=V, root_type=root_type)
         call eos%pressure(&
            n=n, V=V, T=T, P=P_dummy, dPdV=dPdV, dPdT=dPdT, dPdn=dPdn &
            )
         call eos%entropy_residual_vt(&
            n=n, V=V, T=T, Sr=Sr_v, SrV=SrV_v, SrT=SrT_v, Srn=Srn_v &
            )

         dVdP = 1 / dPdV
         dVdT = -dPdT / dPdV
         dVdn = -dPdn / dPdV
      else
         call eos%volume(n=n, P=P, T=T, V=V, root_type=root_type)
         call eos%entropy_residual_vt(n=n, V=V, T=T, Sr=Sr_v)
      end if

      nt = sum(n)
      Z = P * V / (nt * R * T)

      if (present(Sr)) Sr = Sr_v + nt * R * log(Z)

      if (present(SrP)) SrP = SrV_v * dVdP + nt*R*(1.0_pr / P + dVdP / V)

      if (present(SrT)) SrT = SrT_v + SrV_v * dVdT + nt*R*(dVdT / V - 1.0_pr/T)

      if (present(Srn)) Srn = &
         Srn_v + SrV_v * dVdn + R*log(Z) + nt*R*(dVdn / V - 1.0_pr/nt)
   end subroutine entropy_residual_pt

   subroutine internal_energy_residual_pt(&
      eos, n, P, T, root_type, Ur, UrP, UrT, Urn &
      )
      !! Calculate residual internal energy given pressure and temperature.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! eos = PengRobinson76(Tc, Pc, w)
      !!
      !! n = [1.0_pr, 1.0_pr]
      !! T = 300.0_pr
      !! P = 1.0_pr
      !!
      !! call eos%internal_energy_residual_pt(&
      !!    n, P, T, root_type="stable", Ur=Ur, UrP=UrP, UrT=UrT, Urn=Urn &
      !!    )
      !! ```
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]
      character(len=*), intent(in) :: root_type
      !! Desired root-type to solve. Options are:
      !! `["liquid", "vapor", "stable"]`
      real(pr), optional, intent(out) :: Ur !! Residual internal energy [bar L]
      real(pr), optional, intent(out) :: UrP !! \(\frac{dU^r}}{dP}\)
      real(pr), optional, intent(out) :: UrT !! \(\frac{dU^r}}{dT}\)
      real(pr), optional, intent(out) :: Urn(size(n)) !! \(\frac{dU^r}}{dn}\)

      real(pr) :: Ur_v, UrT_v, UrV_v, Urn_v(size(n))
      real(pr) :: V, dPdV, dPdT, dPdn(size(n)), P_dummy
      real(pr) :: dVdP, dVdT, dVdn(size(n))

      logical :: dp, dT, dn, derivs_present

      dp = present(UrP)
      dt = present(UrT)
      dn = present(Urn)

      derivs_present = dp .or. dt .or. dn

      if (derivs_present) then
         call eos%volume(n=n, P=P, T=T, V=V, root_type=root_type)
         call eos%pressure(&
            n=n, V=V, T=T, P=P_dummy, dPdV=dPdV, dPdT=dPdT, dPdn=dPdn &
            )
         call eos%internal_energy_residual_vt(&
            n=n, V=V, T=T, Ur=Ur_v, UrV=UrV_v, UrT=UrT_v, Urn=Urn_v &
            )

         dVdP = 1 / dPdV
         dVdT = -dPdT / dPdV
         dVdn = -dPdn / dPdV
      else
         call eos%volume(n=n, P=P, T=T, V=V, root_type=root_type)
         call eos%internal_energy_residual_vt(n=n, V=V, T=T, Ur=Ur_v)
      end if

      if (present(Ur)) Ur = Ur_v
      if (present(UrP)) UrP = UrV_v * dVdP
      if (present(UrT)) UrT = UrT_v + UrV_v * dVdT
      if (present(Urn)) Urn = Urn_v + UrV_v * dVdn
   end subroutine internal_energy_residual_pt

   subroutine Cv_residual_pt(eos, n, P, T, root_type, Cv)
      !! Calculate residual heat capacity at constant volume given pressure and
      !! temperature.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! eos = PengRobinson76(Tc, Pc, w)
      !!
      !! n = [1.0_pr, 1.0_pr]
      !! T = 300.0_pr
      !! P = 1.0_pr
      !!
      !! call eos%Cv_residual_pt(&
      !!    n, P, T, root_type="stable", Cv=Cv &
      !!    )
      !! ```
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]
      character(len=*), intent(in) :: root_type
      !! Desired root-type to solve. Options are:
      !! `["liquid", "vapor", "stable"]`
      real(pr), intent(out) :: Cv
      !! Residual heat capacity at constant volume [bar L K^-1]

      real(pr) :: V

      call eos%volume(n=n, P=P, T=T, V=V, root_type=root_type)
      call eos%Cv_residual_vt(n=n, V=V, T=T, Cv=Cv)
   end subroutine Cv_residual_pt

   subroutine Cp_residual_pt(eos, n, P, T, root_type, Cp)
      !! Calculate residual heat capacity at constant pressure given pressure
      !! and temperature.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! eos = PengRobinson76(Tc, Pc, w)
      !!
      !! n = [1.0_pr, 1.0_pr]
      !! T = 300.0_pr
      !! P = 1.0_pr
      !!
      !! call eos%Cp_residual_pt(&
      !!    n, P, T, root_type="stable", Cp=Cp &
      !!    )
      !! ```
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]
      character(len=*), intent(in) :: root_type
      !! Desired root-type to solve. Options are:
      !! `["liquid", "vapor", "stable"]`
      real(pr), intent(out) :: Cp
      !! Residual heat capacity at constant pressure [bar L K^-1]

      real(pr) :: V

      call eos%volume(n=n, P=P, T=T, V=V, root_type=root_type)
      call eos%Cp_residual_vt(n=n, V=V, T=T, Cp=Cp)
   end subroutine Cp_residual_pt

   ! ==========================================================================
   ! Excess properties
   ! --------------------------------------------------------------------------
   subroutine ln_activity_coefficient(&
      eos, n, P, T, root_type, lngamma, dlngammadP, dlngammadT, dlngammadn &
      )
      !! Calculate natural logarithm of activity coefficients and its
      !! derivatives given pressure and temperature.
      !!
      !! # Examples
      !!
      !! ```fortran ! eos = PengRobinson76(Tc, Pc, w)
      !!
      !! n = [1.0_pr, 1.0_pr] ! T = 300.0_pr ! P = 1.0_pr
      !!
      !! call eos%ln_activity_coefficient(&
      !!    n, P, T, root_type="stable", &
      !!    lngamma=lngamma, dlngammadP=dlngammadP, &
      !!    dlngammadT=dlngammadT, dlngammadn=dlngammadn &
      !!    )
      !! ```
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]
      character(len=*), intent(in) :: root_type
      !! Desired root-type to solve. Options are:
      !! `["liquid", "vapor", "stable"]`
      real(pr), intent(out), optional :: lngamma(size(n))
      !! Natural logarithm of activity coefficient
      real(pr), intent(out), optional :: dlngammadP(size(n))
      !! \(\frac{d\ln\gamma}}{dP}\)
      real(pr), intent(out), optional :: dlngammadT(size(n))
      !! \(\frac{d\ln\gamma}}{dT}\)
      real(pr), intent(out), optional :: dlngammadn(size(n),size(n))
      !! \(\frac{d\ln\gamma}}{dn}\)

      ! Mixture properties
      real(pr) :: lnPhi(size(n)), dlnPhidT(size(n))
      real(pr) :: dlnPhidn(size(n),size(n))
      real(pr) :: dPdV, dPdn(size(n)), dVdn(size(n))

      ! Pure properties
      real(pr) :: vi(size(n)), vi_temp, npure(size(n))
      real(pr) :: lnPhi_i(size(n)), dlnPhi_i_dT(size(n))
      real(pr) :: lnPhi_i_temp(size(n)), dlnPhi_i_dT_temp(size(n))

      integer :: i

      logical :: gam, dp, dt, dn, present_derivs

      gam = present(lngamma)
      dp = present(dlngammadP)
      dt = present(dlngammadT)
      dn = present(dlngammadn)
      present_derivs = dp .or. dt .or. dn


      ! Call to mixture fugacity coefficient at PT
      ! TODO: maybe more efficient later?
      if (.not. present_derivs) then
         call eos%lnphi_pt(n=n, P=P, T=T, root_type=root_type, lnPhi=lnPhi)
      else
         call eos%lnphi_pt(&
            n=n, P=P, T=T, root_type=root_type, &
            lnPhi=lnPhi, dlnPhidT=dlnPhidT, dlnPhidn=dlnPhidn, &
            dPdV=dPdV, dPdn=dPdn &
            )

         dVdn = -dPdn / dPdV
      end if

      ! Pure components calls
      vi = 0.0_pr
      lnPhi_i = 0.0_pr

      if (dt) then
         dlnPhi_i_dT = 0.0_pr
      end if

      do i=1, size(n)
         npure = 0.0_pr
         npure(i) = 1.0_pr

         if (.not. dt) then
            call eos%lnphi_pt(&
               n=npure, P=P, T=T, V=vi_temp, root_type="stable", &
               lnPhi=lnPhi_i_temp &
               )
         else
            call eos%lnphi_pt(&
               n=npure, P=P, T=T, V=vi_temp, root_type="stable", &
               lnPhi=lnPhi_i_temp, dlnPhidT=dlnPhi_i_dT_temp &
               )
         end if

         vi(i) = vi_temp
         lnPhi_i(i) = lnPhi_i_temp(i)

         if (dt) then
            dlnPhi_i_dT(i) = dlnPhi_i_dT_temp(i)
         end if
      end do

      ! returns
      if (gam) lngamma = lnPhi - lnPhi_i
      if (dt) dlngammadT = dlnPhidT - dlnPhi_i_dT
      if (dp) dlngammadP = (dVdn - vi) / (R * T)
      if (dn) dlngammadn = dlnPhidn
   end subroutine ln_activity_coefficient

   subroutine gibbs_excess(eos, n, P, T, root_type, Ge, GeP, GeT, Gen)
      !! Calculate excess Gibbs energy and its derivatives given pressure and
      !! temperature.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! eos = PengRobinson76(Tc, Pc, w)
      !!
      !! n = [1.0_pr, 1.0_pr]
      !! T = 300.0_pr
      !! P = 1.0_pr
      !!
      !! call eos%gibbs_excess(&
      !!    n, P, T, root_type="stable", &
      !!    Ge=Ge, GeP=GeP, GeT=GeT, Gen=Gen &
      !!    )
      !! ```
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]
      character(len=*), intent(in) :: root_type
      !! Desired root-type to solve. Options are:
      !! `["liquid", "vapor", "stable"]`
      real(pr), intent(out), optional :: Ge
      !! Excess Gibbs energy [bar L]
      real(pr), intent(out), optional :: GeP
      !! \(\frac{dG^E}}{dP}\)
      real(pr), intent(out), optional :: GeT
      !! \(\frac{dG^E}}{dT}\)
      real(pr), intent(out), optional :: Gen(size(n))
      !! \(\frac{dG^E}}{dn}\)

      real(pr) :: lngamma(size(n)), dlngammadP(size(n)), dlngammadT(size(n))
      real(pr) :: dlngammadn(size(n),size(n))

      integer :: j

      logical :: dp, dt, dn, present_derivs

      dp = present(GeP)
      dt = present(GeT)
      dn = present(Gen)
      present_derivs = dp .or. dt .or. dn

      if (.not. present_derivs) then
         call eos%ln_activity_coefficient(&
            n=n, P=P, T=T, root_type=root_type, lngamma=lngamma &
            )
      else
         call eos%ln_activity_coefficient(&
            n=n, P=P, T=T, root_type=root_type, &
            lngamma=lngamma, dlngammadP=dlngammadP, &
            dlngammadT=dlngammadT, dlngammadn=dlngammadn &
            )
      end if

      if (present(Ge)) Ge = R * T * sum(n * lngamma)
      if (present(GeP)) GeP = R * T * sum(n * dlngammadP)
      if (present(GeT)) GeT = R * sum(n * lngamma) + R * T * sum(n * dlngammadT)

      if (present(Gen)) then
         do j=1, size(n)
            Gen(j) = R * T * sum(n * dlngammadn(:,j)) + R * T * lngamma(j)
         end do
      end if
   end subroutine gibbs_excess

   subroutine enthalpy_excess(eos, n, P, T, root_type, He)
      !! Calculate excess enthalpy given pressure and temperature.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! eos = PengRobinson76(Tc, Pc, w)
      !!
      !! n = [1.0_pr, 1.0_pr]
      !! T = 300.0_pr
      !! P = 1.0_pr
      !!
      !! call eos%enthalpy_excess(n, P, T, root_type="stable", He=He)
      !! ```
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]
      character(len=*), intent(in) :: root_type
      !! Desired root-type to solve. Options are:
      !! `["liquid", "vapor", "stable"]`
      real(pr), intent(out) :: He
      !! Excess enthalpy [bar L]

      real(pr) :: dlngammadT(size(n))

      call eos%ln_activity_coefficient(&
         n=n, P=P, T=T, root_type=root_type, dlngammadT=dlngammadT &
         )

      He = -R * T**2 * sum(n * dlngammadT)
   end subroutine enthalpy_excess

   subroutine volume_excess(eos, n, P, T, root_type, Ve)
      !! Calculate excess volume given pressure and temperature.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! eos = PengRobinson76(Tc, Pc, w)
      !!
      !! n = [1.0_pr, 1.0_pr]
      !! T = 300.0_pr
      !! P = 1.0_pr
      !!
      !! call eos%volume_excess(n, P, T, root_type="stable", Ve=Ve)
      !! ```
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]
      character(len=*), intent(in) :: root_type
      !! Desired root-type to solve. Options are:
      !! `["liquid", "vapor", "stable"]`
      real(pr), intent(out) :: Ve !! Excess volume [L]

      real(pr) :: dlngammadP(size(n))

      call eos%ln_activity_coefficient(&
         n=n, P=P, T=T, root_type=root_type, dlngammadP=dlngammadP &
         )

      Ve = R * T * sum(n * dlngammadP)
   end subroutine volume_excess

   subroutine entropy_excess(eos, n, P, T, root_type, Se)
      !! Calculate excess entropy given pressure and temperature.
      !!
      !! # Examples
      !!
      !! ```fortran ! eos = PengRobinson76(Tc, Pc, w)
      !!
      !! n = [1.0_pr, 1.0_pr] ! T = 300.0_pr ! P = 1.0_pr
      !!
      !! call eos%entropy_excess(n, P, T, root_type="stable", Se=Se)
      !! ```
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]
      character(len=*), intent(in) :: root_type
      !! Desired root-type to solve. Options are:
      !! `["liquid", "vapor", "stable"]`
      real(pr), intent(out) :: Se
      !! Excess entropy [bar L K^-1]

      real(pr) :: lngamma(size(n)), dlngammadT(size(n))

      call eos%ln_activity_coefficient(&
         n=n, P=P, T=T, root_type=root_type, &
         lngamma=lngamma, dlngammadT=dlngammadT &
         )

      Se = -R * T * sum(n * dlngammadT) - R * sum(n * lngamma)
   end subroutine entropy_excess

   subroutine helmholtz_excess(eos, n, P, T, root_type, Ae)
      !! Calculate excess Helmholtz energy given pressure and temperature.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! eos = PengRobinson76(Tc, Pc, w)
      !!
      !! n = [1.0_pr, 1.0_pr]
      !! T = 300.0_pr
      !! P = 1.0_pr
      !!
      !! call eos%helmholtz_excess(n, P, T, root_type="stable", Ae=Ae)
      !! ```
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]
      character(len=*), intent(in) :: root_type
      !! Desired root-type to solve. Options are:
      !! `["liquid", "vapor", "stable"]`
      real(pr), intent(out) :: Ae
      !! Excess Helmholtz energy [bar L]

      real(pr) :: lngamma(size(n))

      call eos%ln_activity_coefficient(&
         n=n, P=P, T=T, root_type=root_type, lngamma=lngamma &
         )

      Ae = R * T * sum(n * lngamma)
   end subroutine helmholtz_excess
   ! ==========================================================================
   ! Excess properties
   ! --------------------------------------------------------------------------
   real(pr) function Psat_pure(eos, ncomp, T)
      !! Calculation of saturation pressure of a pure component using the
      !! secant method.
      class(ArModel), intent(in) :: eos !! Model that will be used
      integer, intent(in) :: ncomp
      !! Number of component in the mixture from which the saturation pressure
      !! will be calculated
      real(pr), intent(in) :: T !! Temperature [K]

      real(pr) :: P1, P2
      real(pr) :: f1, f2

      real(pr) :: n(size(eos))

      n = 0
      n(ncomp) = 1

      P1 = 0.5
      P2 = 1

      do while(abs(diff(P2)) > 1e-5)
         f1 = diff(P1)
         f2 = diff(P2)
         Psat_pure = (P1 * f2 - P2 * f1)/(f2 - f1)
         P1 = P2
         P2 = Psat_pure
      end do
   contains
      real(pr) function diff(P)
         real(pr), intent(in) :: P
         real(pr) :: V_l, V_v
         real(pr) :: phi_v(size(eos)), phi_l(size(eos))
         call eos%lnphi_pt(n, P=P, T=T, V=V_v, lnPhi=phi_v, root_type="vapor")
         call eos%lnphi_pt(n, P=P, T=T, V=V_l, lnPhi=phi_l, root_type="liquid")
         diff = phi_v(ncomp) - phi_l(ncomp)
      end function diff
   end function Psat_pure
end module yaeos__models_ar
