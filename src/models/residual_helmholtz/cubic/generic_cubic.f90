module yaeos__models_ar_genericcubic
   use yaeos__constants, only: pr
   use yaeos__models_ar, only: ArModel
   use yaeos__substance, only: Substances
   implicit none

   type, abstract :: AlphaFunction
      !! Abstract derived type that describe the required
      !! procedure for an alpha function.
   contains
      procedure(abs_alpha), deferred :: alpha
   end type

   type, abstract :: CubicMixRule
      !! Abstract derived type that describe the required
      !! procedure for a mixing rule on a Cubic EoS
   contains
      procedure(abs_Dmix), deferred :: Dmix
      procedure(abs_Bmix), deferred :: Bmix
      procedure(abs_D1mix), deferred :: D1mix
   end type

   type, extends(ArModel) :: CubicEoS
      !! # Cubic Equation of State.
      !!
      !! Generic Cubic Equation of State as defined by Michelsen and Mollerup
      !! with a \(\delta_1\) parameter that is not constant, 
      !! and a \(\delta_2\) parameter that depends on it. In the case of a 
      !! two parameter EoS like PengRobinson the \(\delta_1\) is the same for
      !! all components so it can be considered as a constant instead of a 
      !! variable. The expression of the Equation is:
      !!
      !! \[
      !!   P = \frac{RT}{V-B} 
      !!       - \frac{D(T_r)}{(V+B\Delta_1)(V+B\Delta_2)}
      !! \]
      class(CubicMixRule), allocatable :: mixrule
         !! # CubicMixRule derived type.
         !! Uses the abstract derived type `CubicMixRule` to define the
         !! mixing rule that the CubicEoS will use. It includes internally
         !! three methods to calculate the corresponding parameters for the
         !! Cubic EoS: `Dmix`, `Bmix` and `D1mix`.
         !! 
         !! # Examples
         !! ## Calculation of the B parameter.
         !! ```fortran
         !! use yaeos, only: CubicEoS, PengRobinson76
         !! type(CubicEoS) :: eos
         !! eos = PengRobinson76(tc, pc, w)
         !! call eos%mixrule%Bmix(n, eos%b, B, dBi, dBij)
         !! ```
         !! ## Calculation of the D parameter.
         !! ```fortran
         !! use yaeos, only: CubicEoS, PengRobinson76
         !! type(CubicEoS) :: eos
         !! eos = PengRobinson76(tc, pc, w)
         !!
         !! ! The mixing rule takes the `a` parameters of the components so
         !! ! they should be calculated externally
         !! call eos%alpha%alpha(Tr, a, dadt, dadt2)
         !! a = a * eos%ac
         !! dadt = dadt * eos%ac / eos%components%Tc
         !! dadt = dadt * eos%ac / eos%components%Tc**2
         !! ! Calculate parameter
         !! call eos%mixrule%Dmix(n, T, a, dadt, dadt2, D, dDdT, dDdT2, dDi, dDidT, dDij)
         !! ```
         !! ## Calculation of the D1 parameter.
         !! ```fortran
         !! use yaeos, only: CubicEoS, PengRobinson76
         !! type(CubicEoS) :: eos
         !! eos = PengRobinson76(tc, pc, w)
         !! call eos%mixrule%D1mix(n, eos%del1, D1, dD1i, dD1ij)
         !! ```
      class(AlphaFunction), allocatable :: alpha 
         !! # AlphaFunction derived type.
         !! Uses the abstract derived type `AlphaFunction` to define the
         !! Alpha function that the CubicEoS will use. The Alpha function
         !! receives the reduced temperature and returns the values of alpha
         !! and its derivatives, named `a`, `dadt` and `dadt2` respectively.
         !!
         !! # Examples
         !! ## Callign the AlphaFunction of a setted up model.
         !! ```fortran
         !! use yaeos, only: CubicEoS, PengRobinson76
         !! 
         !! type(CubicEoS) :: eos
         !! eos = PengRobinson76(tc, pc, w)
         !! call eos%alpha%alpha(Tr, a, dadt, dadt2)
         !! ```
      real(pr), allocatable :: ac(:) !! Attractive critical parameter
      real(pr), allocatable :: b(:) !! Repulsive parameter
      real(pr), allocatable :: del1(:) !! \(\delta_1\) paramter
      real(pr), allocatable :: del2(:) !! \(\delta_2\) paramter
   contains
      procedure :: residual_helmholtz => GenericCubic_Ar
      procedure :: get_v0 => v0
      procedure :: volume => volume
   end type

   abstract interface
      subroutine abs_alpha(self, Tr, a, dadt, dadt2)
         import AlphaFunction, pr
         class(AlphaFunction), intent(in) :: self
         real(pr), intent(in) :: Tr(:)
         real(pr), intent(out) :: a(:), dadt(:), dadt2(:)
      end subroutine

      subroutine abs_Dmix(self, n, T, &
         ai, daidt, daidt2, &
         D, dDdT, dDdT2, dDi, dDidT, dDij&
         )
         import CubicMixRule, pr
         class(CubicMixRule), intent(in) :: self
         real(pr), intent(in) :: T, n(:)
         real(pr), intent(in) :: ai(:), daidt(:), daidt2(:)
         real(pr), intent(out) :: D, dDdT, dDdT2, dDi(:), dDidT(:), dDij(:, :)
      end subroutine
      
      subroutine abs_Bmix(self, n, bi, B, dBi, dBij)
         import CubicMixRule, pr
         class(CubicMixRule), intent(in) :: self
         real(pr), intent(in) :: n(:)
         real(pr), intent(in) :: bi(:)
         real(pr), intent(out) :: B, dBi(:), dBij(:, :)
      end subroutine
      subroutine abs_D1mix(self, n, d1i, D1, dD1i, dD1ij)
         import pr, CubicMixRule
         class(CubicMixRule), intent(in) :: self
         real(pr), intent(in) :: n(:)
         real(pr), intent(in) :: d1i(:)
         real(pr), intent(out) :: D1
         real(pr), intent(out) :: dD1i(:)
         real(pr), intent(out) :: dD1ij(:, :)
      end subroutine abs_D1mix
   end interface

contains

   subroutine GenericCubic_Ar(&
      self, n, V, T, Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2&
      )
      !! Residual Helmholtz Energy for a generic Cubic Equation of State.
      !!
      !! Calculates the residual Helmholtz Energy for a generic Cubic EoS as
      !! defined by Michelsen and Møllerup:
      !!
      !! \[
      !!   P = \frac{RT}{V-b} 
      !!       - \frac{a_c\alpha(T_r)}{(V+b\delta_1)(V+b\delta_2)}
      !! \]
      !!
      !! This routine assumes that the \(\delta_1\) is not a constant parameter
      !! (as it uses to be in classical Cubic EoS) to be compatible with the
      !! three parameter EoS RKPR where \(delta_1\) is not a constant and
      !! has its own mixing rule.
      !!
      use yaeos__constants, only: R
      class(CubicEoS), intent(in) :: self
      real(pr), intent(in) :: n(:) !! Number of moles
      real(pr), intent(in) :: v !! Volume [L]
      real(pr), intent(in) :: t !! Temperature [K]

      real(pr), optional, intent(out) :: ar !! Residual Helmholtz
      real(pr), optional, intent(out) :: arv !! \(\frac{dAr}{dV}\)
      real(pr), optional, intent(out) :: ArT !! \(\frac{dAr}{dT}\)
      real(pr), optional, intent(out) :: artv !! \(\frac{d^2Ar}{dTdV}\)
      real(pr), optional, intent(out) :: arv2 !! \(\frac{d^2Ar}{dV^2}\)
      real(pr), optional, intent(out) :: ArT2 !! \(\frac{d^2Ar}{dT^2}\)
      real(pr), optional, intent(out) :: Arn(size(n)) !! \(\frac{dAr}{dn_i}\)
      real(pr), optional, intent(out) :: ArVn(size(n)) !! \(\frac{d^2Ar}{dVdn_i}\)
      real(pr), optional, intent(out) :: ArTn(size(n)) !! \(\frac{d^2Ar}{dTdn_i}\)
      real(pr), optional, intent(out) :: Arn2(size(n), size(n)) !! \(\frac{d^2Ar}{dn_{ij}}\)


      real(pr) :: Bmix, dBi(size(n)), dBij(size(n), size(n))
      real(pr) :: D, dDi(size(n)), dDij(size(n), size(n)), dDidT(size(n)), dDdT, dDdT2

      real(pr) :: totn
      real(pr) d1, dD1i(size(n)), dD1ij(size(n), size(n))
      real(pr) :: auxD2, fD1, fBD1, fVD1, fD1D1
      real(pr) d2

      real(pr) :: f, g, fv, fB, gv, fv2, gv2, AUX, FFB, FFBV, FFBB

      real(pr) :: Tr(size(n)), a(size(n)), dadt(size(n)), dadt2(size(n))


      integer :: i, j, nc

      nc = size(n)
      TOTN = sum(n)

     
      Tr = T/self%components%Tc
     
      ! ========================================================================
      ! Attractive parameter and derivatives
      ! ------------------------------------------------------------------------
      call self%alpha%alpha(Tr, a, dadt, dadt2)
      a = self%ac * a
      dadt = self%ac * dadt / self%components%Tc
      dadt2 = self%ac * dadt2 / self%components%Tc**2
      
      ! ========================================================================
      ! Mixing rules
      ! ------------------------------------------------------------------------
      call self%mixrule%D1mix(n, self%del1, D1, dD1i, dD1ij)
      call self%mixrule%Bmix(n, self%b, Bmix, dBi, dBij)
      call self%mixrule%Dmix(&
         n, T, a, dadt, dadt2, D, dDdT, dDdT2, dDi, dDidT, dDij&
         )
      D2 = (1._pr - D1)/(1._pr + D1)

      ! ========================================================================
      ! Main functions defined by Møllerup
      ! The f's and g's used here are for Ar, not F (reduced Ar)
      ! This requires to multiply by R all g, f
      ! ------------------------------------------------------------------------
      f = log((V + D1*Bmix)/(V + D2*Bmix))/Bmix/(D1 - D2)
      g = R*log(1 - Bmix/V)
      fv = -1/((V + D1*Bmix)*(V + D2*Bmix))
      fB = -(f + V*fv)/Bmix
      gv = R*Bmix/(V*(V - Bmix))
      fv2 = (-1/(V + D1*Bmix)**2 + 1/(V + D2*Bmix)**2)/Bmix/(D1 - D2)
      gv2 = R*(1/V**2 - 1/(V - Bmix)**2)

      ! DERIVATIVES OF f WITH RESPECT TO DELTA1
      auxD2 = (1 + 2/(1 + D1)**2)
      fD1 = (1/(V + D1*Bmix) + 2/(V + D2*Bmix)/(1 + D1)**2) - f*auxD2
      fD1 = fD1/(D1 - D2)
      fBD1 = -(fB*auxD2 + D1/(V + D1*Bmix)**2 + 2*D2/(V + D2*Bmix)**2/(1 + D1)**2)
      fBD1 = fBD1/(D1 - D2)
      fVD1 = -(fV*auxD2 + 1/(V + D1*Bmix)**2 + 2/(V + D2*Bmix)**2/(1 + D1)**2)/(D1 - D2)
      fD1D1 = 4*(f - 1/(V + D2*Bmix))/(1 + D1)**3 + Bmix*(-1/(V + D1*Bmix)**2 &
            + 4/(V + D2*Bmix)**2/(1 + D1)**4) - 2*fD1*(1 + 2/(1 + D1)**2)
            fD1D1 = fD1D1/(D1 - D2)

      AUX = R*T/(V - Bmix)
      FFB = TOTN*AUX - D*fB
      FFBV = -TOTN*AUX/(V - Bmix) + D*(2*fv + V*fv2)/Bmix
      FFBB = TOTN*AUX/(V - Bmix) - D*(2*f + 4*V*fv + V**2*fv2)/Bmix**2

      ! ========================================================================
      ! Reduced Helmholtz Energy and derivatives
      ! ------------------------------------------------------------------------
      if (present(Ar)) Ar = -TOTN*g*T - D*f
      if (present(ArV)) ArV = -TOTN*gv*T - D*fv
      if (present(ArV2)) ArV2 = -TOTN*gv2*T - D*fv2

      if (present(Arn))  Arn(:)  = -g*T + FFB*dBi(:) - f*dDi(:) - D*fD1 * dD1i(:)
      if (present(ArVn)) ArVn(:) = -gv*T + FFBV*dBi(:) - fv*dDi(:) - D*fVD1*dD1i(:)
      if (present(ArTn)) ArTn(:) = -g + (TOTN*AUX/T - dDdT*fB)*dBi(:) - f*dDidT(:) - dDdT*fD1*dD1i(:)

      if (present(Arn2)) then
         do i = 1, nc
            do j = 1, i
               Arn2(i, j) = AUX*(dBi(i) + dBi(j)) - fB*(dBi(i)*dDi(j) + dBi(j)*dDi(i)) &
                  + FFB*dBij(i, j) + FFBB*dBi(i)*dBi(j) - f*dDij(i, j)
               Arn2(i, j) = Arn2(i, j) - D*fBD1*(dBi(i)*dD1i(j) + dBi(j)*dD1i(i)) &
                        - fD1*(dDi(i)*dD1i(j) + dDi(j)*dD1i(i)) &
                        - D*fD1*dD1ij(i, j) - D*fD1D1*dD1i(i)*dD1i(j)
               Arn2(j, i) = Arn2(i, j)
            end do
         end do
      end if

      ! TEMPERATURE DERIVATIVES
      if (present(ArT))  ArT = -TOTN*g - dDdT*f
      if (present(ArTV)) ArTV = -TOTN*gv - dDdT*fV
      if (present(ArT2)) ArT2 = -dDdT2*f
   end subroutine GenericCubic_Ar

   function v0(self, n, p, t)
      !! Cubic EoS volume initializer.
      !! For a Cubic Equation of State, the covolume calculated with the mixing
      !! rule is a good estimate for the initial volume solver on the liquid
      !! region.
      class(CubicEoS), intent(in) :: self
      real(pr), intent(in) :: n(:), p, t
      real(pr) :: v0

      real(pr) :: dbi(size(n)), dbij(size(n), size(n))
      call self%mixrule%Bmix(n, self%b, v0, dbi, dbij)
   end function

   subroutine volume(eos, n, P, T, V, root_type)
      !! # Cubic EoS volume solver
      !! Volume solver optimized for Cubic Equations of State.
      !! 
      !! @warn
      !! This routine intends to use the analyitical solution of the cubic
      !! equation, but due to errors in the solutions it is not used. And
      !! the general volume solver by Michelsen is used instead.
      !! @endwarn
      !!
      !! # Description
      !! Cubic equations can be analytically solved. Using an anallytical
      !! solution provides the best possible solution in terms of speed and
      !! precision. This subroutine uses the modified cardano method proposed
      !! by Rosendo.
      !!
      !! # Examples
      !!
      !! ```fortran
      !!  use yaeos, only: CubicEoS, PengRobinson
      !!  type(CubicEoS) :: eos
      !!
      !!  eos = PengRobinson(tc, pc, w)
      !!  ! Possible roots to solve
      !!  call eos%volume(n, P, T, V, "liquid")
      !!  call eos%volume(n, P, T, V, "vapor")
      !!  call eos%volume(n, P, T, V, "stable")
      !! ```
      !!
      !! # References
      !! 
      !! - [1] "Thermodynamic Models: Fundamental and Computational Aspects",
      !!  Michael L. Michelsen, Jørgen M. Mollerup. 
      !!  Tie-Line Publications, Denmark (2004)
      !! [doi](http://dx.doi.org/10.1016/j.fluid.2005.11.032)
      !!
      !! - [2] "A Note on the Analytical Solution of Cubic Equations of State 
      !! in Process Simulation", Rosendo Monroy-Loperena 
      !! [doi](https://dx.doi.org/10.1021/ie2023004)
      use yaeos__constants, only: R
      use yaeos__math_linalg, only: cubic_roots, cubic_roots_rosendo
      use yaeos__models_solvers, only: volume_michelsen
      class(CubicEoS), intent(in) :: eos
      real(pr), intent(in) :: n(:), P, T
      real(pr), intent(out) :: V
      character(len=*), intent(in) :: root_type

      real(pr) :: z(size(n))
      real(pr) :: cp(4), rr(3)
      complex(pr) :: cr(3)
      integer :: flag

      real(pr) :: V_liq, V_vap
      real(pr) :: Ar, AT_Liq, AT_Vap

      real(pr) :: Bmix, dBi(size(n)), dBij(size(n), size(n))
      real(pr) :: D, dDi(size(n)), dDij(size(n), size(n)), dDidT(size(n)), dDdT, dDdT2
      real(pr) :: D1, D2, dD1i(size(n)), dD1ij(size(n), size(n))
      real(pr) :: Tr(size(n))
      real(pr) :: a(size(n)), dadt(size(n)), dadt2(size(n))
      real(pr) :: totn

      call volume_michelsen(eos, n, P, T, V, root_type)
      return

      totn = sum(n)
      z = n/totn
      Tr = T/eos%components%Tc
      ! ========================================================================
      ! Attractive parameter and derivatives
      ! ------------------------------------------------------------------------
      call eos%alpha%alpha(Tr, a, dadt, dadt2)
      a = eos%ac * a
      dadt = eos%ac * dadt / eos%components%Tc
      dadt2 = eos%ac * dadt2 / eos%components%Tc**2
      
      ! ========================================================================
      ! Mixing rules
      ! ------------------------------------------------------------------------
      call eos%mixrule%D1mix(z, eos%del1, D1, dD1i, dD1ij)
      call eos%mixrule%Bmix(z, eos%b, Bmix, dBi, dBij)
      call eos%mixrule%Dmix(&
         z, T, a, dadt, dadt2, D, dDdT, dDdT2, dDi, dDidT, dDij&
         )
      D2 = (1._pr - D1)/(1._pr + D1)

      cp(1) = -P
      cp(2) = -P*Bmix*(D1 + D2 - 1) + R*T
      cp(3) = -P*(D1*D2*Bmix**2 - D1*Bmix**2 - D2*Bmix**2) + R*T*Bmix*(D1+D2) - D
      cp(4) = P*D1*D2*Bmix**3 + R*T *D1*D2*Bmix**2 + D*Bmix

      ! call cubic_roots(cp, rr, cr, flag)
      ! call cubic_roots_rosendo(cp, rr, cr, flag)

      select case(flag)
         case(-1)
            V_liq = rr(1)
            V_vap = rr(3)
            if (V_liq < 0) V_liq = V_vap
         case(1)
            V_liq = rr(1)
            V_vap = rr(1)
      end select

      select case(root_type)
         case("liquid")
            V = V_liq
         case("vapor")
            V = V_vap
         case("stable")
            ! AT is something close to Gr(P,T)
            call eos%residual_helmholtz(z, V_liq, T, Ar=Ar)
            AT_Liq = (Ar + V_liq*P)/(T*R) - sum(z)*log(V_liq)

            call eos%residual_helmholtz(z, V_vap, T, Ar=Ar)
            AT_Vap = (Ar + V_vap*P)/(T*R) - sum(z)*log(V_vap)

            if (AT_liq <= AT_vap) then
               V = V_liq
            else
               V = V_vap
            end if
      end select

      V = totn * V 
   end subroutine
end module
