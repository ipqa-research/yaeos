module generic_cubic
   !! Implementation of the generic cubic Equation of State.
   !!
   !! This module implements a generic cubic based on the definition given
   !! by Michelsen and Møllerup. 
   
   use constants, only: pr, R
   use hyperdual_mod
   use yaeos_interfaces, only: pures_property, abs_cubic_mix
   use ar_models, only: set_ar_function
   use yaeos_thermo_properties, only: vinit

   implicit none

   private

   public :: set_functions, a_res

   procedure(pures_property), pointer :: attractive_parameter !! Attractive parameter function
   procedure(pures_property), pointer :: repulsive_parameter !! Repulsive parameter functoin
   procedure(pures_property), pointer :: delta1_parameter !! \(\delta_1\) parameter
   procedure(pures_property), pointer :: delta2_parameter !! \(\delta_2\) parameter
   procedure(pures_property), pointer :: volume_traslation !! Volume traslation parameter
   procedure(abs_cubic_mix), pointer :: mix !! Mixing function

contains

   subroutine set_functions(a, b, c, d1, d2, mixfun)
      !! Setter of Generic Cubic EoS prcoedures.
      !!
      !! This subroutine receives the desired subroutines to be
      !! used in the calculation of each parameter, and assings the 
      !! corresponding pointers to them. As well as the mixing rule
      !! subroutine to be used and assigns volume initializer function
      procedure(pures_property) :: a
      procedure(pures_property) :: b
      procedure(pures_property) :: c
      procedure(pures_property) :: d1
      procedure(pures_property) :: d2
      procedure(abs_cubic_mix) :: mixfun

      attractive_parameter => a
      repulsive_parameter => b
      delta1_parameter => d1
      delta2_parameter => d2
      volume_traslation => c
      mix => mixfun

      vinit => v0
   end subroutine
      
   pure subroutine a_res(z, v, t, ar)
      !! Generic Cubic EoS Ar function
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
      ! TODO: Update it to also mix delta 1 and 2
      call mix( &
               z, v, t, &
               a_pures, b_pures, c_pures, &
               a, b, c &
      )

      del1 = del1_pures(1)
      del2 = del2_pures(1)

      ar = (&
            -sum(z)*log(1.0_pr - b/v) - a/(R*t*b)*1.0_pr/(del1 - del2) &
            *log((1.0_pr + del1*b/v)/(1.0_pr + del2*b/v))&
      ) ! * R * t
   end subroutine

   pure function v0(z, p, t)
      !! Volume initializer function.
      !!
      !! The classic cubic equation uses the mixture's covolume as
      !! an initializer
      real(pr), intent(in) :: z(:)
      real(pr), intent(in) :: p
      real(pr), intent(in) :: t
      real(pr) :: v0

      type(hyperdual) :: z_d(size(z)), p_d, t_d, b(size(z))

      z_d = z
      p_d = p
      t_d = t

      call repulsive_parameter(z_d, p_d, t_d, b)
      v0 = sum(z*b%f0)*1.2_pr
   end function
end module

module cubic_eos
   !! Set of suboroutines that are part of the classic Cubic EoS, like
   !! PR and SRK alpha function and constant rest of parameters.
   use constants, only: pr, R
   use hyperdual_mod

   implicit none

   ! Model parameters
   real(pr), allocatable :: ac(:) !! Critical atractive parameter
   real(pr), allocatable :: b(:) !! Repulsive parameter
   real(pr), allocatable :: c(:) !! Volume traslation
   real(pr), allocatable :: k(:) !! alpha k-parameter

   real(pr), allocatable :: del1(:) !! \(\delta_1\) Parameter
   real(pr), allocatable :: del2(:)!! \(\delta_2\) Parameter

   ! Critical constants
   real(pr), allocatable :: pc(:) !! Critical Pressure
   real(pr), allocatable :: tc(:) !! Critical Temperature
   real(pr), allocatable :: w(:) !! Acentric factor

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
      !! Atractive parameter
      !!
      !! We assume as classic the PR/SRK equation for the alpha function:
      !! VdW EoS could be a intermediate way 
      !! \[a = a_c (1 + k (1 - \sqrt{\frac{T}{T_c}})))^2\]
      type(hyperdual), intent(in) :: z(:), v, t
      type(hyperdual), intent(out) :: a(size(z))

      a = ac*(1.0_pr + k*(1.0_pr - sqrt(t/tc)))**2
   end subroutine

   pure subroutine b_classic(z, v, t, b_out)
      !! Repulsive parameter
      !!
      !! The repulsive parameter is held constant.
      type(hyperdual), intent(in) :: z(:), v, t
      type(hyperdual), intent(out) :: b_out(size(z))
      b_out = b
   end subroutine
   
   pure subroutine c_classic(z, v, t, c_out)
      !! Volume traslation
      !!
      !! Most of volume traslation methods keep the parameter constant
      type(hyperdual), intent(in) :: z(:), v, t
      type(hyperdual), intent(out) :: c_out(size(z))
      c_out = c
   end subroutine

   pure subroutine del1_classic(z, v, t, del1_out)
      !! \(\delta_1\) Parameter.
      !!
      !! In most cubic equation of state systems this parameter is
      !! held constant.
      type(hyperdual), intent(in) :: z(:), v, t
      type(hyperdual), intent(out) :: del1_out(size(z))
      del1_out = del1
   end subroutine

   pure subroutine del2_classic(z, v, t, del2_out)
      !! \(\delta_1\) Parameter.
      !!
      !! In most cubic equation of state systems this parameter is
      !! held constant.
      type(hyperdual), intent(in) :: z(:), v, t
      type(hyperdual), intent(out) :: del2_out(size(z))
      del2_out = del2
   end subroutine

end module