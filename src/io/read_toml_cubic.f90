! Commented out sections are because they depend on the new API that's on
! another branch

module toml_cubic
   use constants
   use thermo_io
   use cubic_eos, only: CubicEOS

   implicit none

   interface alloc
      !! Allocate the arrays based on the number of components
      module procedure :: cubic_system_alloc
   end interface alloc

   interface read_system
      !! Read the system from a toml file
      module procedure :: read_cubic_system
   end interface read_system

   interface setup_from_toml
      !! Setup a model from a toml file
      ! module procedure :: setup_from_toml_PengRobinson76
      ! module procedure :: setup_from_toml_PengRobinson78
      ! module procedure :: setup_from_toml_SoaveRedlichKwong
   end interface

contains

   subroutine cubic_system_alloc(self, n)
      class(CubicEOS) :: self
      integer, intent(in) :: n

      allocate (self%tc(n))
      allocate (self%pc(n))
      allocate (self%w(n))

      allocate (self%ac(n))
      allocate (self%b(n))
      allocate (self%k(n))
   end subroutine cubic_system_alloc

   subroutine read_cubic_system(conf_file, system)
      ! Relevant data variables
      character(len=*), intent(in) :: conf_file
      class(CubicEOS), intent(in out) :: system

      character(len=:), allocatable :: thermo_model
      character(len=:), allocatable :: name
      character(len=:), allocatable :: component_file

      integer :: i, n_components

      call read_system(conf_file, system%ArModel)

      thermo_model = system%thermo_model
      n_components = size(system%z)

      call alloc(system, n_components)

      do i = 1, n_components
         name = system%names(i)

         ! Look for component inside the main input file, if it's not there
         ! open the database file
         call get_value(main_table, name, component_table, requested=.false.)
         if (.not. associated(component_table)) then
            component_file = trim(database_path)//trim(name)//".toml"
            print *, component_file
            call toml_load(temp_table, component_file)
            call get_value(temp_table, name, component_table)
         end if

         ! Reading pure parameters
         call get_value(component_table, "pure", params_table)

         call component_table%get_keys(keys)

         call get_value(params_table, "Tc", system%tc(i))
         call get_value(params_table, "Pc", system%pc(i))
         call get_value(params_table, "w", system%w(i))

         ! Reading model parameters
         if (use_parameters) then
            call get_value(component_table, thermo_model, params_table)

            call get_value(params_table, "a", system%ac(i))
            call get_value(params_table, "b", system%b(i))
            call get_value(params_table, "k", system%k(i))
         end if
      end do
      ! =======================================================================
   end subroutine read_cubic_system

!    subroutine setup_from_toml_PengRobinson76(toml_file, system)
!       character(len=*), intent(in) :: toml_file
!       type(PengRobinson76) :: system
! 
!       call read_system(toml_file, system%ArModel)
!       call read_system(toml_file, system%CubicEOS)
! 
!       if (.not. use_parameters) then
!          call setup(system, from_constants=.true.)
!       end if
!    end subroutine setup_from_toml_PengRobinson76
! 
!    subroutine setup_from_toml_PengRobinson78(toml_file, system)
!       character(len=*), intent(in) :: toml_file
!       type(PengRobinson78) :: system
! 
!       call read_system(toml_file, system%ArModel)
!       call read_system(toml_file, system%CubicEOS)
!    end subroutine setup_from_toml_PengRobinson78
! 
!    subroutine setup_from_toml_SoaveRedlichKwong(toml_file, system)
!       character(len=*), intent(in) :: toml_file
!       type(SoaveRedlichKwong) :: system
! 
!       call read_system(toml_file, system%ArModel)
!       call read_system(toml_file, system%CubicEOS)
!    end subroutine setup_from_toml_SoaveRedlichKwong

end module toml_cubic
