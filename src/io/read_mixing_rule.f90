module toml_mixingrules
   use constants, only: pr, database_path
   use thermo_io
   use models, only : ArModel
   use mixing_rules, only: ClassicVdW

   implicit none

   interface alloc
      !! Allocate the arrays based on the number of components
      module procedure :: alloc_classicvdw
   end interface alloc

   interface read_system
      !! Read the system from a toml file
      module procedure :: read_classicvdw
   end interface read_system

contains

    subroutine alloc_classicvdw(self, n)
        class(ClassicVdW), intent(in out) :: self
        integer, intent(in) :: n

        allocate(self%kij(n, n))
        allocate(self%lij(n, n))
    end subroutine alloc_classicvdw

    subroutine read_classicvdw(conf_file, system, mix_rule)
        ! Relevant data variables
        character(len=*), intent(in) :: conf_file
        class(ArModel), intent(in out) :: system
        class(ClassicVdW), intent(in out) :: mix_rule

        integer :: i, j

        associate(&
            model => system%thermo_model, n_components => size(system), &
            names => system%names &
        )
            call alloc(mix_rule, n_components)

            do i=1, n_components
                do j=i,n_components
                end do
            end do
               
        end associate
    end subroutine
end module
