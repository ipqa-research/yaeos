module alpha_mathias_copeman
   use yaeos_models_ar_cubic_alphas, only: AlphaFunction
   use yaeos_constants, only: pr
   implicit none

   type, extends(AlphaFunction) :: MathiasCopeman
      real(pr), allocatable :: c1(:)
      real(pr), allocatable :: c2(:)
      real(pr), allocatable :: c3(:)
   contains
      procedure :: alpha => alpha
   end type

contains

   subroutine alpha(self, Tr, a, dadt, dadt2)
      class(MathiasCopeman), intent(in) :: self
      real(pr), intent(in) :: Tr(:)
      real(pr), intent(out) :: a(:), dadt(:), dadt2(:)

      real(pr) :: sqrt_Tr(size(Tr))

      sqrt_Tr = 1 - sqrt(Tr)

      associate(c1 => self%c1, c2 => self%c2, c3 => self%c3)
         a = (1 + c1 * (sqrt_Tr) + c2 * (sqrt_Tr) + c3 * (sqrt_Tr))**2
         dadt = (c1 + c2 + c3) * (c1*(sqrt(Tr) - 1) &
                 + c2*(sqrt(Tr) - 1) + c3*(sqrt(Tr) - 1) - 1)/sqrt(Tr)
         dadt2 = (1.0_pr/2.0_pr) * (&
                 c1**2 + 2*c1*c2 + 2*c1*c3 &
                 + c1 + c2**2 + 2*c2*c3 +  c2 + c3**2 + c3)/Tr**(3.0_pr/2.0_pr)
      end associate
   end subroutine
end module

module new_alpha_example
    use yaeos, only: pr, PengRobinson76, CubicEoS
    use alpha_mathias_copeman, only: MathiasCopeman
    implicit none
contains
    subroutine main
        use example_models, only:binary_PR76
        use yaeos, only: fugacity_vt
        class(CubicEoS), allocatable :: eos
        type(MathiasCopeman) :: alpha

        real(pr) :: n(2), v, t
        real(pr) :: lnfug(2)

        ! Get the example PR76 binary model
        eos = binary_PR76()

        ! Define the new alpha function parameters
        alpha%c1 = [0.5, 0.3]
        alpha%c2 = [0.6, 0.8]
        alpha%c3 = [0.6, 0.8]

        ! Replace the original alpha
        eos%alpha = alpha

        n = [0.3, 0.7]
        v = 1
        t = 150
        call fugacity_vt(eos, n, v, t, lnfug=lnfug)

        print *, lnfug
    end subroutine
end module
