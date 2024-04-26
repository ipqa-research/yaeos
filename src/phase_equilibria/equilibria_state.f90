module yaeos_equilibria_equilibria_state
   use yaeos_constants, only: pr
   implicit none

   type :: EquilibriaState
      integer :: iters !! Iterations needed to reach the state
      real(pr), allocatable :: y(:) !! Light-phase molar fractions
      real(pr), allocatable :: x(:) !! Heavy-phase molar fractions
      real(pr) :: Vx !! Liquid volume [L/mol]
      real(pr) :: Vy !! Vapor volume [L/mol]
      real(pr) :: t !! Temperature [K]
      real(pr) :: p !! Pressure [bar]
      real(pr) :: beta !! Mole fraction of light-phase
   end type
end module
