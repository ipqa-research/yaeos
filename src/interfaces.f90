module yaeos_interfaces
   use hyperdual_mod

   implicit none
   private
   public :: dual_property
   public :: abs_cubic_mix
   public :: abs_cubic_prop_mix
   public :: pures_property

   abstract interface
     pure subroutine dual_property(z, v, t, property)
        import hyperdual
        type(hyperdual), intent(in) :: z(:), v, t
        type(hyperdual), intent(out) :: property
     end subroutine dual_property
   end interface
   
   abstract interface
     pure subroutine pures_property(z, v, t, property)
        import hyperdual
        type(hyperdual), intent(in) :: z(:), v, t
        type(hyperdual), intent(out) :: property(size(z))
     end subroutine pures_property
   end interface

   abstract interface
     pure subroutine abs_cubic_mix(z, v, t, a, b, c, amix, bmix, cmix)
         import hyperdual
         type(hyperdual), intent(in) :: z(:), v, t, a(size(z)), b(size(z)), c(size(z))
         type(hyperdual), intent(out) :: amix, bmix, cmix
     end subroutine
   end interface

   abstract interface
     pure subroutine abs_cubic_prop_mix(z, v, t, prop, prop_mix)
         import hyperdual
         type(hyperdual), intent(in) :: z(:), v, t, prop(size(z))
         type(hyperdual), intent(out) :: prop_mix
     end subroutine
   end interface
end module