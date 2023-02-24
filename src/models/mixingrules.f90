module mixing_rules
   !!  Mixing rules module
   use constants, only: pr
   use hyperdual_mod

   implicit none

   ! ===========================================================================
   ! Mixing Rules Derived Types
   ! ---------------------------------------------------------------------------
   type, abstract :: CubicMixingRule
      !! Cubic EoS Mixing Rule
      !! Contains a single subroutine as an attribute that receives the
      !! set of parameters of pure compounds and returns the mixture's
      !! parameters and their derivatives
   contains
      procedure(abs_cubic_mix), deferred :: mix
   end type

   type, extends(CubicMixingRule) :: ClassicVdW
      !! # Classic Mixing rules
      !!
      !! ## attractive term:
      !! \(a_{mix} = \sum_i z_i \sum_j z_j \sqrt{a_i a_j} (1-k_{ij})\)
      !!
      !! ## Repulsive term:
      !! \(b_{mix} = \sum_i z_i \sum_j z_j \frac{b_i b_j}{2} (1 - l_{ij})\)
      !!
      !! ## Volume shift:
      !! \(c_{mix} = \sum_i z_i c_i \)
      !!
      !! ## \(\delta_1 \) mix:
      !! First component value (assumed all \(\delta_1\) equal)
      !!
      !! ## \(\delta_2\) mix:
      !! First component value (assumed all \(\delta_2\) equal)
      real(pr), allocatable :: kij(:, :) !! \(k_{ij}\) matrix
      real(pr), allocatable :: lij(:, :) !! \(l_{ij}\)
   contains
      procedure :: mix => mix_ClassicVdW
   end type
   ! ===========================================================================

   ! ===========================================================================
   ! Interfaces
   ! ---------------------------------------------------------------------------
   abstract interface
      pure subroutine abs_cubic_mix(mixrule, z, v, t, a, b, c, amix, bmix, cmix)
         import CubicMixingRule
         import hyperdual
         class(CubicMixingRule), intent(in) :: mixrule
         type(hyperdual), intent(in) :: z(:), v, t, a(size(z)), b(size(z)), c(size(z))
         type(hyperdual), intent(out) :: amix, bmix, cmix
      end subroutine

      pure subroutine abs_cubic_prop_mix(mixrule, z, v, t, prop, prop_mix)
         import CubicMixingRule
         import hyperdual
         class(CubicMixingRule), intent(in) :: mixrule
         type(hyperdual), intent(in) :: z(:), v, t, prop(size(z))
         type(hyperdual), intent(out) :: prop_mix
      end subroutine
   end interface

   interface setup_MixingRule
      module procedure :: setup_ClassicVdW
   end interface

   interface mix
      module procedure :: mix_ClassicVdW
   end interface
   ! ===========================================================================

contains

   ! ===========================================================================
   ! Constructors
   ! ---------------------------------------------------------------------------
   subroutine setup_ClassicVdW(mixrule, kij, lij)
      type(ClassicVdW) :: mixrule
      real(pr) :: kij(:, :)
      real(pr) :: lij(:, :)

      mixrule%kij = kij
      mixrule%lij = lij
   end subroutine
   ! ===========================================================================

   ! ===========================================================================
   ! Mixing subroutines
   ! ---------------------------------------------------------------------------
   pure subroutine mix_ClassicVdW(mixrule, z, v, t, a, b, c, amix, bmix, cmix)
      class(ClassicVdW), intent(in) :: mixrule
      type(hyperdual), intent(in) :: z(:), v, t, a(size(z)), b(size(z)), c(size(z))
      type(hyperdual), intent(out) :: amix, bmix, cmix
      
      call a_mix(mixrule, z, v, t, a, amix)
      call b_mix(mixrule, z, v, t, b, bmix)
      call c_mix(mixrule, z, v, t, c, cmix)
   end subroutine

   pure subroutine a_mix_classic(mixrule, z, v, t, a, amix)
      class(ClassicVdW), intent(in) :: mixrule
      type(hyperdual), intent(in) :: z(:), v, t, a(size(z))
      type(hyperdual), intent(out) :: amix

      type(hyperdual) :: aij(size(z), size(z))

      integer :: i, j

      aij = 0.0_pr
      associate (kij => mixrule%kij)
         amix = 0.0_pr
         do j = 1, size(z)
            aij(:, j) = sqrt(a(:) * a(j)) * (1.0_pr - kij(:, j))
            amix = amix + sum(z(:) * z(j) * aij(:, j))
         end do
      end associate
   end subroutine

   pure subroutine b_mix_classic(mixrule, z, v, t, b, bmix)
      class(ClassicVdW), intent(in) :: mixrule
      type(hyperdual), intent(in) :: z(:), v, t, b(size(z))
      type(hyperdual), intent(out) :: bmix

      type(hyperdual) :: bij(size(z), size(z))

      integer :: i

      associate (lij => mixrule%lij)
         bmix = 0.0_pr
         do i = 1, size(z)
            bij(:, i) = (b(:) + b(i))/2.0_pr*(1.0_pr - lij(:, i))
            bmix = bmix + sum(z(:)*z(i)*bij(:, i))
         end do

         bmix = bmix/sum(z)
      end associate
   end subroutine

   pure subroutine c_mix_classic(mixrule, z, v, t, c, cmix)
      class(ClassicVdW), intent(in) :: mixrule
      type(hyperdual), intent(in) :: z(:), v, t, c(size(z))
      type(hyperdual), intent(out) :: cmix

      cmix = sum(z*c)
   end subroutine
   ! ===========================================================================

end module mixing_rules
