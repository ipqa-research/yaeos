module generic_cubic
   use constants, only: pr, R
   use hyperdual_mod
   use yaeos_interfaces, only: pures_property, abs_cubic_mix
   use ar_models, only: set_ar_function

   procedure(pures_property), pointer :: attractive_parameter
   procedure(pures_property), pointer :: repulsive_parameter
   procedure(pures_property), pointer :: delta1_parameter
   procedure(pures_property), pointer :: delta2_parameter
   procedure(pures_property), pointer :: volume_traslation
   procedure(abs_cubic_mix), pointer :: mix

contains

   subroutine set_parameters_functions(a, b, d1, d2, c)
      procedure(pures_property) :: a
      procedure(pures_property) :: b
      procedure(pures_property) :: d1
      procedure(pures_property) :: d2
      procedure(pures_property) :: c
      attractive_parameter => a
      repulsive_parameter => b
      delta1_parameter => d1
      delta2_parameter => d2
      volume_traslation => c
   end subroutine
      
   subroutine set_mixrule(mixfun)
      procedure(abs_cubic_mix) :: mixfun
      mix => mixfun
   end subroutine

   pure subroutine a_res(z, v, t, ar)
      ! Generic Cubic EoS 
      type(hyperdual), intent(in) :: z(:), v, t
      type(hyperdual), intent(out) :: ar

      type(hyperdual), dimension(size(z)) :: a_pures, b_pures, c_pures
      type(hyperdual), dimension(size(z)) ::  del1_pures, del2_pures
      type(hyperdual) :: a, b, c, del1, del2

      ! Pure components parameters
      call attractive_parameter(z, v, t, a_pures)
      call repulsive_parameter(z, v, t, b_pures)

      call delta1_parameter(z, v, t, del1_pures)
      call delta2_parameter(z, v, t, del2_pures)

      ! Call the model's mixing rule
      call mix( &
               z, v, t, &
               a_pures, b_pures, c_pures, &
               a, b, c &
      )

      ar = (&
            -sum(z)*log(1.0_pr - b/v) - a/(R*t*b)*1.0_pr/(del1 - del2) &
            *log((1.0_pr + del1*b/v)/(1.0_pr + del2*b/v))&
      )! * R * t
   end subroutine
end module

module cubic_eos
   use constants, only: pr, R
   use hyperdual_mod

   implicit none

   ! Model parameters
   real(pr), allocatable :: ac(:)
   real(pr), allocatable :: b(:)
   real(pr), allocatable :: c(:)
   real(pr), allocatable :: k(:)

   real(pr), allocatable :: del1(:)
   real(pr), allocatable :: del2(:)

   ! Critical constants
   real(pr), allocatable :: pc(:)
   real(pr), allocatable :: tc(:)
   real(pr), allocatable :: w(:)

contains

   ! Allocators
   subroutine alloc(n)
      integer, intent(in) :: n
      allocate(ac(n), b(n), c(n), k(n), del1(n), del2(n), pc(n), tc(n), w(n))
   end subroutine
   
   subroutine destroy()
      deallocate(ac, b, c, k, del1, del2, pc, tc, w)
   end subroutine

   pure subroutine a_classic(z, v, t, a)
      type(hyperdual), intent(in) :: z(:), v, t
      type(hyperdual), intent(out) :: a(size(z))

      a = ac*(1.0_pr + k*(1.0_pr - sqrt(t/tc)))**2
   end subroutine

   pure subroutine b_classic(z, v, t, b)
      !! Repulsive parameter
      !! The repulsive parameter is held constant.
      type(hyperdual), intent(in) :: z(:), v, t
      type(hyperdual), intent(out) :: b(size(z))
   end subroutine
   
   pure subroutine c_classic(z, v, t, c)
      type(hyperdual), intent(in) :: z(:), v, t
      type(hyperdual), intent(out) :: c(size(z))
   end subroutine

   pure subroutine del1_classic(z, v, t, del1)
      type(hyperdual), intent(in) :: z(:), v, t
      type(hyperdual), intent(out) :: del1(size(z))
   end subroutine

   pure subroutine del2_classic(z, v, t, del2)
      type(hyperdual), intent(in) :: z(:), v, t
      type(hyperdual), intent(out) :: del2(size(z))
   end subroutine

end module
