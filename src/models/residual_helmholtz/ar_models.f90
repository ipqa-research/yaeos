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
      use yaeos__constants, only: pr, R
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
      use yaeos__constants, only: R
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
      use yaeos__constants, only: R
      class(ArModel), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: V !! Volume [L]
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(out) :: Cp !! heat capacity P constant [bar L / K]

      real(pr) :: Cv, P, dPdT, dPdV, n_t

      n_t = sum(n)

      call Cv_residual_vt(eos, n, V, T, Cv)
      call pressure(eos, n, V, T, P, dPdV=dPdV, dPdT=dPdT)

      Cp = Cv - T * dPdT**2 / dPdV - n_t*R
   end subroutine Cp_residual_vt

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
