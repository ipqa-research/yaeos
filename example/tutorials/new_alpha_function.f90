! In this example we will define a new alpha function for a cubic Equation of
! state. In this case it is the Mathias-Copeman alpha function.
module alpha_mathias_copeman
   ! The base AlphaFunction structure is defined in the
   ! yaeos_models_ar_cubic_alphas module
   use yaeos_models_ar_cubic_alphas, only: AlphaFunction
   use yaeos_constants, only: pr
   implicit none

   ! We will extend that structure to include the parameters c1, c2 and c3
   ! used in the function.
   type, extends(AlphaFunction) :: MathiasCopeman
      real(pr), allocatable :: c1(:)
      real(pr), allocatable :: c2(:)
      real(pr), allocatable :: c3(:)
   contains
      ! And define its alpha functionality with the subroutine defined in the
      ! module.
      procedure :: alpha => alpha
   end type MathiasCopeman

contains

   subroutine alpha(self, Tr, a, dadt, dadt2)
      !! MathiasCopeman alpha function definition
      class(MathiasCopeman), intent(in) :: self
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
   end subroutine alpha
end module alpha_mathias_copeman

program new_alpha_example
   use yaeos__example_tools, only: methane_butane_pr76
   use yaeos, only: fugacity_vt, pressure
   use yaeos, only: pr, PengRobinson76, CubicEoS, QMR
   use alpha_mathias_copeman, only: MathiasCopeman
   type(CubicEoS) :: eos
   type(QMR) :: mr
   type(MathiasCopeman) :: alpha

   real(pr) :: n(2), v, t, P
   integer :: i

   ! Get the example PR76 binary model
   eos = methane_butane_pr76()

   ! Define the new alpha function parameters
   alpha%c1 = [0.49258, 0.84209]
   alpha%c2 = [0.0, -0.46406]
   alpha%c3 = [0.0, 0.84619]
   
   n = [0.3, 0.7]
   v = 2
   t = 150

   call pressure(eos, n, V, T, P=P)
   print *, "Peng-Robinson76:", P
   
   ! Replace the original alpha
   deallocate(eos%alpha) ! Remove the already defined alpha
   eos%alpha = alpha     ! assign the new defined alpha

   call pressure(eos, n, V, T, P=P)
   print *, "Peng-Robinson76-MC:", P


end program new_alpha_example
