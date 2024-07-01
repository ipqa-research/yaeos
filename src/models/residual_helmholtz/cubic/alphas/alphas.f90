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
   end type

   type, extends(AlphaFunction) :: AlphaRKPR
      !! RKPR \(\alpha\) function
      !! \[
      !! \alpha(T_r) = \left(\frac{3}{2 + T_r}\right)^k
      !! \]
      real(pr), allocatable :: k(:) !! \(k\) parameter.
   contains
      procedure :: alpha => alpha_rkpr
   end type

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
   end subroutine

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
   end subroutine
end module