module yaeos__models_ar_genericcubic
   !! # Generic cubic equation of state
   !!
   !! # Description
   !! Base types and the generic implementation of a cubic equation of state.
   !! The module defines the abstract [[AlphaFunction]] and [[CubicMixRule]]
   !! types, and the [[CubicEoS]] model that combines them. Specific cubic
   !! equations are provided in
   !! [[yaeos__models_ar_cubic_implementations(module)]].
   !! Defines the alpha-function and mixing-rule contracts, evaluates generic
   !! cubic Helmholtz expressions, and provides model configuration methods.
   !!
   !! # Examples
   !! ```fortran
   !! model = PengRobinson76(Tc, Pc, omega)
   !! call model%set_alpha(alpha_function)
   !! ```
   use yaeos__constants, only: pr
   use yaeos__models_ar, only: ArModel
   use yaeos__substance, only: Substances
   implicit none

   type, abstract :: AlphaFunction
      !! # AlphaFunction
      !! Abstract contract for temperature corrections to cubic attraction.
      !!
      !! # Description
      !! Implementations return component-wise \(\alpha(T_r)\) and its first
      !! two derivatives with respect to reduced temperature. `CubicEoS`
      !! converts those derivatives to absolute-temperature derivatives.
      !!
      !! # Examples
      !! See `AlphaSoave`, `AlphaRKPR`, and `AlphaMathiasCopeman`.
   contains
      procedure(abs_alpha), deferred :: alpha
   end type AlphaFunction

   type, abstract :: CubicMixRule
      !! # CubicMixRule
      !! Abstract contract for cubic-EOS mixture parameter calculations.
      !!
      !! # Description
      !! Concrete mixing rules implement attractive, repulsive, and cubic
      !! denominator parameter mixing, including derivatives needed for
      !! residual properties and phase equilibrium.
      !!
      !! # Examples
      !! See `QMR`, `CMR`, `MHV`, and `HV_NRTL` for implementations.
      logical :: is_D_ddlc = .false.
      !! Mixing rule D parameter depends on volume (density)
      logical :: dn2 = .false.
      !! Calculate second-order mole-number derivatives
   contains
      procedure(abs_Dmix), deferred :: Dmix
      procedure(abs_Bmix), deferred :: Bmix
      procedure(abs_D1mix), deferred :: D1mix
   end type CubicMixRule

   type, extends(ArModel) :: CubicEoS
      !! # CubicEoS
      !! Generic cubic equation of state with configurable alpha and mixing rules.
      !!
      !! # Description
      !! Generic Cubic Equation of State as defined by Michelsen and Mollerup
      !! with a \(\delta_1\) parameter that is not constant,
      !! and a \(\delta_2\) parameter that depends on it. In the case of a
      !! two parameter EoS like PengRobinson the \(\delta_1\) is the same for
      !! all components so it can be considered as a constant instead of a
      !! variable. The expression of the Equation is:
      !!
      !! \[
      !!   P = \frac{RT}{V-B}
      !!       - \frac{D(T_r)}{(V+B\delta_1)(V+B\delta_2)}
      !! \]
      !!
      !! # Examples
      !! ```fortran
      !! model = PengRobinson76(Tc, Pc, omega)
      !! call model%pressure(n, V, T, P)
      !! ```
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
      !! dadt2 = dadt2 * eos%ac / eos%components%Tc**2
      !! ! Calculate parameter
      !! call eos%mixrule%Dmix(&
      !!    n, V, T, a, dadt, dadt2, D, &
      !!    dDdV, dDdT, dDdV2, dDdT2, dDi, dDdTV, dDidV, dDidT, dDij &
      !!    )
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
      !! ## Calling the AlphaFunction of a configured model.
      !! ```fortran
      !! use yaeos, only: CubicEoS, PengRobinson76
      !!
      !! type(CubicEoS) :: eos
      !! eos = PengRobinson76(tc, pc, w)
      !! call eos%alpha%alpha(Tr, a, dadt, dadt2)
      !! ```
      real(pr), allocatable :: ac(:) !! Pure-component attractive parameters at \(T_c\), \(a_{c,i}\) [bar L^2/mol^2]
      real(pr), allocatable :: b(:) !! Pure-component covolumes \(b_i\) [L/mol]
      real(pr), allocatable :: del1(:) !! \(\delta_1\) parameter
      real(pr), allocatable :: del2(:) !! \(\delta_2\) parameter
   contains
      procedure :: residual_helmholtz => GenericCubic_Ar
      procedure :: get_v0 => v0
      procedure :: volume => volume
      procedure :: set_delta1 => set_delta1
      procedure :: set_mixrule => set_mixrule
      procedure :: set_alpha => set_alpha
      procedure :: calculate_attractive_parameters => calculate_attractive_parameters
   end type CubicEoS

   abstract interface
      subroutine abs_alpha(self, Tr, a, dadt, dadt2)
         !! # abs_alpha
         !!
         !! # Description
         !! Evaluate a cubic alpha function and its reduced-temperature derivatives.
         import AlphaFunction, pr
         class(AlphaFunction), intent(in) :: self !! Alpha-function implementation.
         real(pr), intent(in) :: Tr(:) !! Reduced temperatures, one per component.
         real(pr), intent(out) :: a(:) !! Alpha-function values.
         real(pr), intent(out) :: dadt(:) !! First reduced-temperature derivative of alpha.
         real(pr), intent(out) :: dadt2(:) !! Second reduced-temperature derivative of alpha.
      end subroutine abs_alpha

      subroutine abs_Dmix(self, n, V, T, &
         ai, daidt, daidt2, &
         D, &
         dDdV, dDdT, dDdV2, dDdT2, dDi, dDdTV, dDidV, dDidT, dDij &
         )
         !! # abs_Dmix
         !!
         !! # Description
         !! Evaluate the attractive mixture parameter and its state/composition
         !! derivatives for a cubic mixing rule.
         import CubicMixRule, pr
         class(CubicMixRule), intent(in) :: self !! Mixing-rule implementation.
         real(pr), intent(in) :: V !! Total volume [L].
         real(pr), intent(in) :: T !! Temperature [K].
         real(pr), intent(in) :: n(:) !! Component mole amounts [mol].
         real(pr), intent(in) :: ai(:) !! Pure-component attractive parameters.
         real(pr), intent(in) :: daidt(:) !! Temperature derivatives of pure-component attractive parameters.
         real(pr), intent(in) :: daidt2(:) !! Second temperature derivatives of pure-component attractive parameters.
         real(pr), intent(out) :: D !! Attractive mixture parameter.
         real(pr), intent(out) :: dDdV !! Volume derivative of the attractive mixture parameter.
         real(pr), intent(out) :: dDdT !! Temperature derivative of the attractive mixture parameter.
         real(pr), intent(out) :: dDdV2 !! Second volume derivative of the attractive mixture parameter.
         real(pr), intent(out) :: dDdT2 !! Second temperature derivative of the attractive mixture parameter.
         real(pr), intent(out) :: dDdTV !! Mixed temperature/volume derivative of the attractive mixture parameter.
         real(pr), intent(out) :: dDi(:) !! Mole-number gradient of the attractive mixture parameter.
         real(pr), intent(out) :: dDidV(:) !! Volume derivative of the attractive mole-number gradient.
         real(pr), intent(out) :: dDidT(:) !! Temperature derivative of the attractive mole-number gradient.
         real(pr), intent(out) :: dDij(:, :) !! Mole-number Hessian of the attractive mixture parameter.
      end subroutine abs_Dmix

      subroutine abs_Bmix(self, n, bi, B, dBi, dBij)
         !! # abs_Bmix
         !!
         !! # Description
         !! Evaluate the repulsive mixture parameter and its composition derivatives.
         import CubicMixRule, pr
         class(CubicMixRule), intent(in) :: self !! Mixing-rule implementation.
         real(pr), intent(in) :: n(:) !! Component mole amounts [mol].
         real(pr), intent(in) :: bi(:) !! Pure-component covolumes.
         real(pr), intent(out) :: B !! Mixture covolume \(B\) [L].
         real(pr), intent(out) :: dBi(:) !! Mole-number gradient of the mixture covolume.
         real(pr), intent(out) :: dBij(:, :) !! Mole-number Hessian of the mixture covolume.
      end subroutine abs_Bmix

      subroutine abs_D1mix(self, n, d1i, D1, dD1i, dD1ij)
         !! # abs_D1mix
         !!
         !! # Description
         !! Evaluate the first denominator parameter and its composition derivatives.
         import pr, CubicMixRule
         class(CubicMixRule), intent(in) :: self !! Mixing-rule implementation.
         real(pr), intent(in) :: n(:) !! Component mole amounts [mol].
         real(pr), intent(in) :: d1i(:) !! Pure-component \(\delta_1\) values.
         real(pr), intent(out) :: D1 !! Mixture \(\delta_1\) parameter.
         real(pr), intent(out) :: dD1i(:) !! Gradient of `D1` with respect to mole amounts.
         real(pr), intent(out) :: dD1ij(:, :) !! Hessian of `D1` with respect to mole amounts.
      end subroutine abs_D1mix
   end interface

contains

   subroutine GenericCubic_Ar(&
      self, n, V, T, Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2&
      )
      !! # GenericCubic_Ar
      !!
      !! # Description
      !! Residual Helmholtz Energy for a generic Cubic Equation of State.
      !!
      !! Calculates the residual Helmholtz Energy for a generic Cubic EoS as
      !! defined by Michelsen and Møllerup:
      !!
      !! \[
      !!   P = \frac{nRT}{V-B}
      !!       - \frac{D(T,n)}{(V+B\delta_1)(V+B\delta_2)}
      !! \]
      !!
      !! This routine assumes that the \(\delta_1\) is not a constant parameter
      !! (as it uses to be in classical Cubic EoS) to be compatible with the
      !! three parameter EoS RKPR where \(delta_1\) is not a constant and
      !! has its own mixing rule.
      !!
      use yaeos__models_ar_genericcubic_base, only: &
         generic => GenericCubic_Ar, &
         generic_dddlc => GenericCubic_Ar_dddlc
      use yaeos__constants, only: R
      class(CubicEoS), intent(in) :: self !! Cubic model containing parameters and mixing-rule state.
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


      real(pr) :: B
      real(pr) :: dBi(size(n)) !! Mole-number gradient of the mixture covolume.
      real(pr) :: dBij(size(n), size(n)) !! Mole-number Hessian of the mixture covolume.
      real(pr) :: D !! Attractive mixture parameter.
      real(pr) :: dDi(size(n)) !! Mole-number gradient of the attractive mixture parameter.
      real(pr) :: dDij(size(n), size(n)) !! Mole-number Hessian of the attractive mixture parameter.
      real(pr) :: dDidT(size(n)) !! Temperature derivative of the attractive mole-number gradient.
      real(pr) :: dDdT !! Temperature derivative of the attractive mixture parameter.
      real(pr) :: dDdT2 !! Second temperature derivative of the attractive mixture parameter.
      real(pr) :: dDdV !! Volume derivative of the attractive mixture parameter.
      real(pr) :: dDdV2 !! Second volume derivative of the attractive mixture parameter.
      real(pr) :: dDdTV !! Mixed temperature/volume derivative of the attractive mixture parameter.
      real(pr) :: dDidV(size(n)) !! Volume derivative of the attractive mole-number gradient.

      real(pr) :: totn !! Total number of moles.
      real(pr) d1, dD1i(size(n)), dD1ij(size(n), size(n))


      real(pr) :: Tr(size(n))
      real(pr) :: a(size(n))
      real(pr) :: dadt(size(n))
      real(pr) :: dadt2(size(n))


      integer :: nc !! Number of components.

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
      call self%mixrule%Bmix(n, self%b, B, dBi, dBij)
      dDdV = 0
      dDdV2 = 0
      dDdTV = 0
      dDdTV = 0
      call self%mixrule%Dmix(&
         n=n, V=V, T=T, ai=a, daidt=dadt, daidt2=dadt2, &
         D=D, &
         dDdV=dDdV, dDdV2=dDdV2,dDdT=dDdT, dDdT2=dDdT2, &
         dDdTV=dDdTV, dDi=dDi, dDidV=dDidV, dDidT=dDidT, dDij=dDij &
         )

      if (self%mixrule%is_D_ddlc) then
         call generic_dddlc(&
            n=n, V=V, T=T, &
            B=B, &
            dBi=dBi, dBij=dBij, &
            D=D, &
            dDdV=dDdV, dDdV2=dDdV2,dDdT=dDdT, dDdT2=dDdT2, &
            dDdTV=dDdTV, dDi=dDi, dDidV=dDidV, dDidT=dDidT, dDij=dDij, &
            D1=D1, dD1i=dD1i, dD1ij=dD1ij, &
            Ar=Ar, ArV=ArV, ArT=ArT, ArTV=ArTV, ArV2=ArV2, &
            ArT2=ArT2, Arn=Arn, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2&
            )
      else
         call generic(&
            n=n, V=V, T=T, &
            B=B, dBi=dBi, dBij=dBij, &
            D=D, dDi=dDi, dDij=dDij, dDidT=dDidT, dDdT=dDdT, dDdT2=dDdT2, &
            D1=D1, dD1i=dD1i, dD1ij=dD1ij, &
            Ar=Ar, ArV=ArV, ArT=ArT, ArTV=ArTV, ArV2=ArV2, ArT2=ArT2, Arn=Arn, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2&
            )
      end if
   end subroutine GenericCubic_Ar

   subroutine set_delta1(self, delta1)
      !! # set_delta1
      !!
      !! # Description
      !! Set the \(\delta_1\) parameter of the cubic equation.
      !!
      !! The \(\delta_2\) parameter is recalculated from
      !! \(\delta_2 = (1 - \delta_1)/(1 + \delta_1)\).
      class(CubicEoS) :: self
      real(pr), intent(in) :: delta1(:) !! \(\delta_1\) values
      self%del1 = delta1
      self%del2 = (1._pr - delta1)/(1._pr + delta1)
   end subroutine set_delta1

   subroutine set_mixrule(self, mixrule)
      !! # set_mixrule
      !!
      !! # Description
      !! Set the mixing rule used to calculate mixture attractive parameters.
      class(CubicEoS), intent(in out) :: self !! Model whose mixing rule is replaced.
      class(CubicMixRule), intent(in) :: mixrule !! Mixing rule
      if (allocated(self%mixrule)) deallocate(self%mixrule)
      self%mixrule = mixrule
   end subroutine set_mixrule

   subroutine set_alpha(self, alpha)
      !! # set_alpha
      !!
      !! # Description
      !! Set the alpha function used to calculate the attractive parameter.
      class(CubicEoS), intent(in out) :: self !! Model whose alpha function is replaced.
      class(AlphaFunction), intent(in) :: alpha !! Alpha function
      if (allocated(self%alpha)) deallocate(self%alpha)
      self%alpha = alpha
   end subroutine set_alpha

   subroutine calculate_attractive_parameters(self, T, a, dadt, dadt2)
      !! # calculate_attractive_parameters
      !! Pure-component attractive parameters and their temperature derivatives.
      !!
      !! # Description
      !! Evaluates the alpha function at \(T_r = T/T_{c,i}\) and returns
      !! \(a_i = a_{c,i}\alpha_i\), \(da_i/dT\), and \(d^2a_i/dT^2\).
      class(CubicEoS), intent(in out) :: self !! Model
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(out) :: a(:) !! Attractive parameters \(a_i(T) = a_{c,i}\alpha_i\) [bar L^2/mol^2]
      real(pr), intent(out) :: dadt(:) !! \(da_i/dT\) [bar L^2/(mol^2 K)]
      real(pr), intent(out) :: dadt2(:) !! \(d^2a_i/dT^2\) [bar L^2/(mol^2 K^2)]

      real(pr) :: Tr(size(a)) !! Reduced temperatures for each pure component

      Tr = T/self%components%Tc
      call self%alpha%alpha(Tr, a, dadt, dadt2)

      a = self%ac * a
      dadt = self%ac * dadt / self%components%Tc
      dadt2 = self%ac * dadt2 / self%components%Tc**2
   end subroutine calculate_attractive_parameters

   function v0(self, n, p, t)
      !! Cubic EoS volume initializer.
      !! For a Cubic Equation of State, the covolume calculated with the mixing
      !! rule is a good estimate for the initial volume solver on the liquid
      !! region.
      class(CubicEoS), intent(in) :: self !! Cubic model with component covolumes and mixing rule.
      real(pr), intent(in) :: n(:) !! Component mole amounts [mol].
      real(pr), intent(in) :: p !! Pressure [bar].
      real(pr), intent(in) :: t !! Temperature [K].
      real(pr) :: v0

      real(pr) :: dbi(size(n))
      real(pr) :: dbij(size(n), size(n))
      call self%mixrule%Bmix(n, self%b, v0, dbi, dbij)
   end function v0

   subroutine volume(eos, n, P, T, V, root_type)
      !! # Cubic EoS volume solver
      !! Volume solver optimized for Cubic Equations of State.
      !!
      !! # Description
      !! Uses the solver proposed by Michelsen.
      !!
      !! # Examples
      !!
      !! ```fortran
      !!  use yaeos, only: CubicEoS, PengRobinson76
      !!  type(CubicEoS) :: eos
      !!
      !!  eos = PengRobinson76(tc, pc, w)
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
      class(CubicEoS), intent(in) :: eos !! Cubic EOS and mixing-rule parameters.
      real(pr), intent(in) :: n(:) !! Component mole amounts [mol].
      real(pr), intent(in) :: P !! Pressure [bar].
      real(pr), intent(in) :: T !! Temperature [K].
      real(pr), intent(out) :: V !! Selected total volume [L].
      character(len=*), intent(in) :: root_type !! Requested `liquid`, `vapor`, or `stable` root.

      real(pr) :: z(size(n)) !! Feed or overall composition.
      real(pr) :: cp(4)
      real(pr) :: rr(3)
      complex(pr) :: cr(3)
      integer :: flag !! Root/solver classification code.

      real(pr) :: V_liq
      real(pr) :: V_vap
      real(pr) :: Ar
      real(pr) :: AT_Liq
      real(pr) :: AT_Vap

      real(pr) :: Bmix !! Mixture covolume.
      real(pr) :: dBi(size(n)) !! Mole-number gradient of the mixture covolume.
      real(pr) :: dBij(size(n), size(n)) !! Mole-number Hessian of the mixture covolume.
      real(pr) :: D !! Attractive mixture parameter.
      real(pr) :: dDi(size(n)) !! Mole-number gradient of the attractive mixture parameter.
      real(pr) :: dDij(size(n), size(n)) !! Mole-number Hessian of the attractive mixture parameter.
      real(pr) :: dDidT(size(n)) !! Temperature derivative of the attractive mole-number gradient.
      real(pr) :: dDdT !! Temperature derivative of the attractive mixture parameter.
      real(pr) :: dDdT2 !! Second temperature derivative of the attractive mixture parameter.
      real(pr) :: D1
      real(pr) :: D2
      real(pr) :: dD1i(size(n)) !! Mole-number gradient of the cubic denominator parameter.
      real(pr) :: dD1ij(size(n), size(n)) !! Mole-number Hessian of the cubic denominator parameter.
      real(pr) :: Tr(size(n))
      real(pr) :: a(size(n)) !! Alpha-function values.
      real(pr) :: dadt(size(n)) !! First reduced-temperature derivative of alpha.
      real(pr) :: dadt2(size(n)) !! Second reduced-temperature derivative of alpha.
      real(pr) :: totn !! Total number of moles.

      call volume_michelsen(eos, n=n, P=P, T=T, V=V, root_type=root_type)
   end subroutine volume
end module yaeos__models_ar_genericcubic
