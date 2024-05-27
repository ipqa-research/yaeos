module yaeos_equilibria_equilibria_state
   use yaeos_constants, only: pr
   implicit none

   type :: EquilibriaState
      !! Description of a two-phase equilibria state.
      !!
      !! Contains the relevant information of an equilibrium point obtained
      !! from some kind of equilibria calculation.
      character(len=14) :: kind
      !! Kind of point ["bubble", "dew", "liquid-liquid", "split"]
      integer :: iters = 0
      !! Iterations needed to reach the state
      real(pr), allocatable :: y(:)
      !! Light-phase molar fractions
      real(pr), allocatable :: x(:)
      !! Heavy-phase molar fractions
      real(pr) :: Vx
      !! Heavy-phase volume [L/mol]
      real(pr) :: Vy
      !! Light-phase volume [L/mol]
      real(pr) :: T
      !! Temperature [K]
      real(pr) :: P
      !! Pressure [bar]
      real(pr) :: beta
      !! Mole fraction of light-phase
   contains
      private
      procedure, pass :: write => write_EquilibriaState
      generic, public :: write (FORMATTED) => write
   end type EquilibriaState

contains

   subroutine write_EquilibriaState(eq, unit, iotype, v_list, iostat, iomsg)
      class(EquilibriaState), intent(in) :: eq
      integer, intent(in) :: unit
      character(*), optional, intent(in) :: iotype
      integer, optional, intent(in)  :: v_list(:)
      integer, intent(out) :: iostat
      character(*), optional, intent(inout) :: iomsg

      character(*), parameter :: nl = new_line("G")

      write(unit, *) eq%kind, eq%T, eq%P

   end subroutine write_EquilibriaState
end module yaeos_equilibria_equilibria_state
