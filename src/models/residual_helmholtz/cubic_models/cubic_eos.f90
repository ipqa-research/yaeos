module cubic_eos
   use constants, only: pr, R
   use hyperdual_mod
   use ar_models, only: ArModel, dual_property, size, dualderiv, alloc
   use mixing_rules, only: CubicMixingRule

   implicit none

   type, extends(ArModel) :: CubicEOS
      !! Generic Cubic Equation of State
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

      ! MixingRule object
      class(CubicMixingRule), pointer :: mixrule
   end type

   abstract interface
      subroutine dual_cubic_property(model, z, p, v, t, property)
         import hyperdual
         import CubicEOS
         class(CubicEOS) :: model
         type(hyperdual), intent(in) :: z(model%size), p, v, t
         type(hyperdual), intent(out) :: property(model%size)
      end subroutine dual_cubic_property
   end interface

   interface attractive_parameter
      module procedure :: a_classic
   end interface attractive_parameter

   interface repulsive_parameter
      module procedure :: b_classic
   end interface repulsive_parameter

   interface volume_traslation
      module procedure :: c_classic
   end interface volume_traslation

   interface del1_parameter
      module procedure :: del1_classic
   end interface

   interface del2_parameter
      module procedure :: del2_classic
   end interface

   interface ar
      module procedure :: a_res
   end interface

   interface alloc
      module procedure :: alloc_CubicEoS
   end interface alloc
contains
   subroutine alloc_CubicEoS(model, n)
      type(CubicEoS) :: model
      integer :: n

      allocate (model%ac(n))
      allocate (model%b(n))
      allocate (model%c(n))
      allocate (model%k(n))

      allocate (model%del1(n))
      allocate (model%del2(n))

      allocate (model%pc(n))
      allocate (model%tc(n))
      allocate (model%w(n))

      call alloc(model%ArModel, n)
   end subroutine alloc_CubicEoS

   pure subroutine a_classic(model, z, v, t, a)
      type(CubicEOS), intent(in) :: model
      type(hyperdual), intent(in) :: z(size(model)), v, t
      type(hyperdual), intent(out) :: a(size(model))

      associate (ac => model%ac, tc => model%tc, k => model%k)
         a = ac*(1.0_pr + k*(1.0_pr - sqrt(t/tc)))**2
      end associate
   end subroutine

   pure subroutine b_classic(model, z, v, t, b)
      type(CubicEOS), intent(in) :: model
      type(hyperdual), intent(in) :: z(size(model)), v, t

      type(hyperdual), intent(out) :: b(size(model))

      b = model%b
   end subroutine

   pure subroutine c_classic(model, z, v, t, c)
      type(CubicEOS), intent(in) :: model
      type(hyperdual), intent(in) :: z(size(model)), v, t
      type(hyperdual), intent(out) :: c(size(model))

      c = model%c
   end subroutine

   pure subroutine del1_classic(model, z, v, t, del1)
      type(CubicEOS), intent(in) :: model
      type(hyperdual), intent(in) :: z(size(model)), v, t
      type(hyperdual), intent(out) :: del1(size(model))

      del1 = model%del1
   end subroutine

   pure subroutine del2_classic(model, z, v, t, del2)
      type(CubicEOS), intent(in) :: model
      type(hyperdual), intent(in) :: z(size(model)), v, t
      type(hyperdual), intent(out) :: del2(size(model))

      del2 = model%del2
   end subroutine

   pure subroutine a_res(model, z, v, t, ar)
      class(ArModel), intent(in) :: model
      type(hyperdual), intent(in) :: z(size(model)), v, t
      type(hyperdual), intent(out) :: ar

      type(hyperdual), dimension(size(model)) :: a_pures, b_pures, c_pures
      type(hyperdual), dimension(size(model)) ::  del1_pures, del2_pures
      type(hyperdual) :: a, b, c, del1, del2

      select type (model)
      class is (CubicEoS)
         ! Pure components parameters
         call attractive_parameter(model, z, v, t, a_pures)
         call repulsive_parameter(model, z, v, t, b_pures)

         call del1_parameter(model, z, v, t, del1_pures)
         call del2_parameter(model, z, v, t, del2_pures)

         del1 = del1_pures(1)
         del2 = del2_pures(1)

         ! Call the model's mixing rule
         call model%mixrule%mix( &
                  z, v, t, &
                  a_pures, b_pures, c_pures, &
                  a, b, c &
         )

         ar = (&
               -sum(z)*log(1.0_pr - b/v) - a/(R*t*b)*1.0_pr/(del1 - del2) &
               *log((1.0_pr + del1*b/v)/(1.0_pr + del2*b/v))&
         )! * R * t
      end select
   end subroutine
end module
