module yaeos_interfaces
   !-> Set of abstract interfaces that show how the general functions should
   !   behave.
   use hyperdual_mod

   implicit none
   private
   public :: dual_property
   public :: abs_cubic_mix
   public :: abs_cubic_prop_mix
   public :: pures_property
   public :: volume_initalizer

   abstract interface
      pure subroutine dual_property(z, v, t, property)
         !-> Subroutine to calculate a single scalar property, like Ar
         import hyperdual
         type(hyperdual), intent(in) :: z(:), v, t
         type(hyperdual), intent(out) :: property
      end subroutine dual_property

      pure subroutine pures_property(z, v, t, property)
         !-> Subroutine to calculate a mixture property for a set of components,
         ! Like the Cubic EoS attractive parameter
         import hyperdual
         type(hyperdual), intent(in) :: z(:), v, t
         type(hyperdual), intent(out) :: property(size(z))
      end subroutine pures_property

      pure subroutine abs_cubic_mix(z, v, t, a, b, c, amix, bmix, cmix)
         !->  Subroutine to mix a set of parameters for a cubic equation of state
         import hyperdual
         type(hyperdual), intent(in) :: z(:), v, t, a(size(z)), b(size(z)), c(size(z))
         type(hyperdual), intent(out) :: amix, bmix, cmix
      end subroutine

      pure subroutine abs_cubic_prop_mix(z, v, t, prop, prop_mix)
         !->  Subroutine to mix a single property
         import hyperdual
         type(hyperdual), intent(in) :: z(:), v, t, prop(size(z))
         type(hyperdual), intent(out) :: prop_mix
      end subroutine

      pure function volume_initalizer(z, p, t) result(v0)
         !-> Volume initializer routine
         import pr
         real(pr), intent(in) :: z(:)
         real(pr), intent(in) :: p
         real(pr), intent(in) :: t
         real(pr) :: v0
      end function
   end interface
end module
