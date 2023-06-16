module thermo_io
   use constants, only: pr, database_path
   use models, only: ArModel, size
   use tomlf, only: toml_table, toml_array, toml_key, get_value, len
   use tomlf_de, only: toml_load


   implicit none


   ! TOML related variables
   type(toml_table), allocatable :: main_table, temp_table
   type(toml_table), pointer :: system_table, single_composition
   type(toml_table), pointer :: component_table, params_table
   type(toml_array), pointer :: components, composition
   type(toml_key), allocatable :: keys(:)


   ! General system variables
   logical :: use_parameters 
   !! Use model defined parameters instead of getting them from critical 
   !! constants, this variable probably should be inside the main system 
   !! derived type


   interface alloc
      module procedure :: main_system_alloc
   end interface alloc


   interface read_system
      module procedure :: read_main_system
   end interface read_system


contains


   pure function system_size(sys) result(sys_size)
        class(ArModel), intent(in) :: sys
        integer :: sys_size

        sys_size = size(sys%names)
   end function


   subroutine main_system_alloc(self, n)
      type(ArModel) :: self
      integer, intent(in) :: n
      integer :: stat

      allocate(self%z(n), stat=stat)
      allocate(self%names(n), stat=stat)
   end subroutine main_system_alloc


   subroutine read_main_system(conf_file, main_system)
      character(len=*), intent(in) :: conf_file
      type(ArModel), intent(in out) :: main_system

      integer :: i, n, stat
      character(len=:), allocatable :: name
      
      ! ========================================================================
      ! General settings
      call toml_load(main_table, conf_file)
      call get_value(main_table, "system", system_table)
      call get_value(system_table, "use_parameters", use_parameters)

      ! ------------------------------------------------------------------------
      ! Read used model, mixing rule
      call get_value(system_table, "model", main_system%thermo_model)
      call get_value(system_table, "mixing_rule", main_system%mixing_rule)
      
      ! ------------------------------------------------------------------------
      ! Read components
      call get_value( &
         system_table, "components", components, stat=stat, requested=.false. &
      )

      call get_value( &
         system_table, "composition", composition, stat=stat, requested=.false. &
      )

      n = len(composition)
      call alloc(main_system, n)

      do i=1, n
          ! call get_value(components, i, name)
          call get_value(composition, i, single_composition)
          call get_value(single_composition, "name", name)
          call get_value(single_composition, "z", main_system%z(i))

          main_system%names(i) = name
      end do
   end subroutine read_main_system

end module thermo_io