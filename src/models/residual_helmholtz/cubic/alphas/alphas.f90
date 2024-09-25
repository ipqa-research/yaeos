module yaeos__models_ar_cubic_alphas
   !! \(\alpha\) functions defined in the library.
   use yaeos__constants, only: pr
   use yaeos__substance, only: substances
   use yaeos__models_ar_genericcubic, only: CubicEoS, AlphaFunction
   implicit none

   type, extends(AlphaFunction) :: AlphaSoave
      !! Soave \(\alpha\) function.
      !! \( \alpha(T_r) = (1 + k (1 - \sqrt{Tr}))^2 \)
      real(pr), allocatable :: k(:) !! \(k\) parameter.
   contains
      procedure :: alpha !! Alpha function
   end type AlphaSoave

   type, extends(AlphaFunction) :: AlphaRKPR
      !! RKPR \(\alpha\) function
      !! \[
      !! \alpha(T_r) = \left(\frac{3}{2 + T_r}\right)^k
      !! \]
      real(pr), allocatable :: k(:) !! \(k\) parameter.
   contains
      procedure :: alpha => alpha_rkpr
   end type AlphaRKPR


   type, extends(AlphaFunction) :: AlphaMathiasCopeman
      !! Mathias Copeman \(\alpha\) function.
      real(pr), allocatable :: c1(:)
      real(pr), allocatable :: c2(:)
      real(pr), allocatable :: c3(:)
   contains
      procedure :: alpha => alpha_mc
   end type AlphaMathiasCopeman

contains

   subroutine alpha(self, Tr, a, dadt, dadt2)
      !! Soave \(\alpha\) function and it's derivatives.
      class(AlphaSoave), intent(in) :: self
      real(pr), intent(in) :: Tr(:) !! Reduced temperature
      real(pr), intent(out) :: a(:) !! \(\alpha\)
      real(pr), intent(out) :: dadt(:) !! \(\frac{d\alpha}{dT}\)
      real(pr), intent(out) :: dadt2(:)!! \(\frac{d^2\alpha}{dT^2}\)

      associate(k => self%k)
         a = (1 + k*(1 - sqrt(Tr)))**2
         dadT = k*(k*(sqrt(Tr) - 1) - 1)/sqrt(Tr)
         dadT2 = (1.0_pr/2.0_pr)*k*(k + 1)/Tr**(1.5_pr)
      end associate
   end subroutine alpha

   subroutine alpha_rkpr(self, Tr, a, dadt, dadt2)
      class(AlphaRKPR), intent(in) :: self
      real(pr), intent(in) :: Tr(:) !! Reduced temperature
      real(pr), intent(out) :: a(:) !! \(\alpha\)
      real(pr), intent(out) :: dadt(:) !! \(\frac{d\alpha}{dT}\)
      real(pr), intent(out) :: dadt2(:)!! \(\frac{d^2\alpha}{dT^2}\)

      associate(k => self%k)
         a = (3/(2 + Tr))**k
         dadT = -k*a/(2 + Tr)
         dadT2 = -(k + 1)*dadT/(2 + Tr)
      end associate
   end subroutine alpha_rkpr

   subroutine alpha_mc(self, Tr, a, dadt, dadt2)
      !! MathiasCopeman alpha function definition
      class(AlphaMathiasCopeman), intent(in) :: self
      real(pr), intent(in) :: Tr(:)
      real(pr), intent(out) :: a(:), dadt(:), dadt2(:)

      real(pr) :: sqrt_Tr(size(Tr))

      sqrt_Tr = 1 - sqrt(Tr)

      ! The associate statement allows to abreviate the expresions
      associate(c1 => self%c1, c2 => self%c2, c3 => self%c3)
         a = (1 + c1 * (sqrt_Tr) + c2 * (sqrt_Tr) + c3 * (sqrt_Tr))**2
         dadt = (c1 + c2 + c3) * (c1*(sqrt(Tr) - 1) &
            + c2*(sqrt(Tr) - 1) + c3*(sqrt(Tr) - 1) - 1)/sqrt(Tr)
         dadt2 = (1.0_pr/2.0_pr) * (&
            c1**2 + 2*c1*c2 + 2*c1*c3 &
            + c1 + c2**2 + 2*c2*c3 +  c2 + c3**2 + c3)/Tr**(3.0_pr/2.0_pr)
      end associate
   end subroutine alpha_mc

end module yaeos__models_ar_cubic_alphas
