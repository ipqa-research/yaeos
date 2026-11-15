module yaeos__models_ar_cubic_alphas
   !! # Cubic alpha functions
   !! Component alpha-function implementations used by cubic equations of state.
   !!
   !! # Description
   !! The module provides Soave, RKPR, and Mathias-Copeman temperature
   !! corrections. Each type evaluates \(\alpha(T_r)\) and its first two
   !! derivatives with respect to reduced temperature, as required by
   !! `CubicEoS`.
   !!
   !! # Examples
   !! ```fortran
   !! type(AlphaSoave) :: alpha
   !! alpha%k = k_values
   !! call model%set_alpha(alpha)
   !! ```
   use yaeos__constants, only: pr
   use yaeos__substance, only: substances
   use yaeos__models_ar_genericcubic, only: CubicEoS, AlphaFunction
   implicit none

   type, extends(AlphaFunction) :: AlphaSoave
      !! # AlphaSoave
      !! Soave temperature correction for cubic EOS models.
      !!
      !! # Description
      !! Evaluates \(\alpha(T_r)=[1+k(1-\sqrt{T_r})]^2\) componentwise and
      !! supplies its reduced-temperature derivatives.
      !!
      !! # Examples
      !! ```fortran
      !! alpha%k = k_values
      !! call model%set_alpha(alpha)
      !! ```
      real(pr), allocatable :: k(:) !! \(k\) parameter.
   contains
      procedure :: alpha !! Alpha function
   end type AlphaSoave

   type, extends(AlphaFunction) :: AlphaRKPR
      !! # AlphaRKPR
      !! RKPR alpha function for cubic EOS models.
      !!
      !! # Description
      !! Evaluates \(\alpha(T_r)=[3/(2+T_r)]^k\) and its first two
      !! reduced-temperature derivatives componentwise.
      !!
      !! # Examples
      !! ```fortran
      !! alpha%k = k_values
      !! call model%set_alpha(alpha)
      !! ```
      real(pr), allocatable :: k(:) !! \(k\) parameter.
   contains
      procedure :: alpha => alpha_rkpr
   end type AlphaRKPR

   type, extends(AlphaFunction) :: AlphaMathiasCopeman
      !! # AlphaMathiasCopeman
      !! Mathias-Copeman alpha function with component-specific coefficients.
      !!
      !! # Description
      !! For \(T_r \leq 1\):
      !! \[
      !! \alpha = \left[1 + c_1(1-\sqrt{T_r}) + c_2(1-\sqrt{T_r})^2
      !!          + c_3(1-\sqrt{T_r})^3\right]^2
      !! \]
      !! and for \(T_r > 1\) the Soave form with \(c_1\) only,
      !! \(\alpha = [1 + c_1(1-\sqrt{T_r})]^2\).
      !! Each coefficient array has one value per component.
      !!
      !! # Examples
      !! ```fortran
      !! alpha%c1 = c1
      !! alpha%c2 = c2
      !! alpha%c3 = c3
      !! call model%set_alpha(alpha)
      !! ```
      real(pr), allocatable :: c1(:) !! First Mathias-Copeman coefficient by component.
      real(pr), allocatable :: c2(:) !! Second Mathias-Copeman coefficient by component.
      real(pr), allocatable :: c3(:) !! Third Mathias-Copeman coefficient by component.
   contains
      procedure :: alpha => alpha_mc
   end type AlphaMathiasCopeman

contains

   subroutine alpha(self, Tr, a, dadt, dadt2)
      !! Soave \(\alpha\) function and it's derivatives.
      class(AlphaSoave), intent(in) :: self !! Soave alpha parameters by component.
      real(pr), intent(in) :: Tr(:) !! Reduced temperature
      real(pr), intent(out) :: a(:) !! \(\alpha\)
      real(pr), intent(out) :: dadt(:) !! \(\frac{d\alpha}{dT_r}\)
      real(pr), intent(out) :: dadt2(:) !! \(\frac{d^2\alpha}{dT_r^2}\)

      associate(k => self%k)
         a = (1 + k*(1 - sqrt(Tr)))**2
         dadT = k*(k*(sqrt(Tr) - 1) - 1)/sqrt(Tr)
         dadT2 = (1.0_pr/2.0_pr)*k*(k + 1)/Tr**(1.5_pr)
      end associate
   end subroutine alpha

   subroutine alpha_rkpr(self, Tr, a, dadt, dadt2)
      class(AlphaRKPR), intent(in) :: self !! RKPR alpha exponent by component.
      real(pr), intent(in) :: Tr(:) !! Reduced temperature
      real(pr), intent(out) :: a(:) !! \(\alpha\)
      real(pr), intent(out) :: dadt(:) !! \(\frac{d\alpha}{dT_r}\)
      real(pr), intent(out) :: dadt2(:) !! \(\frac{d^2\alpha}{dT_r^2}\)

      associate(k => self%k)
         a = (3/(2 + Tr))**k
         dadT = -k*a/(2 + Tr)
         dadT2 = -(k + 1)*dadT/(2 + Tr)
      end associate
   end subroutine alpha_rkpr

   subroutine alpha_mc(self, Tr, a, dadt, dadt2)
      !! # alpha_mc
      !! Evaluate the Mathias-Copeman alpha function and its reduced-temperature derivatives.
      !!
      !! # Description
      !! Applies the low-temperature polynomial branch for \(T_r\leq1\) and
      !! the high-temperature continuation for \(T_r>1\).
      !!
      !! # Examples
      !! Configure `c1`, `c2`, and `c3` by component and call the type-bound
      !! `alpha` procedure.
      class(AlphaMathiasCopeman), intent(in) :: self !! Component-specific coefficients.
      real(pr), intent(in) :: Tr(:) !! Reduced temperatures.
      real(pr), intent(out) :: a(:) !! Alpha-function values.
      real(pr), intent(out) :: dadt(:) !! First reduced-temperature derivative of alpha.
      real(pr), intent(out) :: dadt2(:) !! Second reduced-temperature derivative of alpha.

      real(pr) :: sqrt_Tr(size(Tr))

      real(pr) :: u(size(Tr))
      real(pr) :: dudt(size(Tr))
      real(pr) :: dudt2(size(Tr))

      sqrt_Tr = 1 - sqrt(Tr)

      ! The associate statement allows to abreviate the expresions
      associate(c1 => self%c1, c2 => self%c2, c3 => self%c3)
         where (Tr > 1)
            a = (1 + c1*(1 - sqrt(Tr)))**2
            dadT = c1*(c1*(sqrt(Tr) - 1) - 1)/sqrt(Tr)
            dadT2 = (1.0_pr/2.0_pr)*c1*(c1 + 1)/Tr**(1.5_pr)
         elsewhere
            u = -c1*(sqrt(Tr) - 1) + c2*(sqrt(Tr) - 1)**2 - c3*(sqrt(Tr) - 1)**3 + 1
            dudt = (1.0_pr/2.0_pr)*(-c1 + 2*c2*(sqrt(Tr) - 1) - 3*c3*(sqrt(Tr) - 1)**2) /sqrt(Tr)
            dudt2 = (1.0_pr/4.0_pr)*(-3*Tr*c3 + c1 + 2*c2 + 3*c3)/Tr**(3.0_pr/2.0_pr)
            
            a = u**2
            dadt = 2*u*dudt
            dadt2 = 2*(dudt**2 + u*dudt2)
         end where
      end associate
   end subroutine alpha_mc

end module yaeos__models_ar_cubic_alphas
