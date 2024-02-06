module yaeos_equilibria_equilibria_state
   use yaeos_constants, only: pr
   implicit none

   type :: EquilibriaState
      integer :: iters !! Iterations needed to reach the state
      real(pr), allocatable :: y(:) !! Vapour molar fractions
      real(pr), allocatable :: x(:) !! Liquid molar fractions
      real(pr) :: t !! Temperature [K]
      real(pr) :: p !! Pressure [bar]
   end type
end module
