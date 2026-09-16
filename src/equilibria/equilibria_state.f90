module yaeos__equilibria_equilibrium_state
   use yaeos__constants, only: pr
   implicit none

   type :: EquilibriumState
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
      procedure, pass :: write => write_EquilibriumState
      generic, public :: write (FORMATTED) => write
   end type EquilibriumState

   type :: MPEquilibriumState
      !! # `MPEquilibriumState`
      !! Type to hold the state of a multiphase equilibrium calculation.
      !!
      !! # Description
      !! This type holds the results of a multiphase equilibrium calculation,
      !! including phase compositions, pressures, and temperatures.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! ```
      !!
      !! # References
      !!
      real(pr), allocatable :: z(:) !! Global composition
      real(pr) :: P !! Pressure
      real(pr) :: T !! Temperature
      integer :: np !! Number of phases
      real(pr), allocatable :: x_l(:,:)  !! Mole fractions of the main phases
      real(pr), allocatable :: w(:)  !! Mole fractions of the reference phase
      real(pr), allocatable :: betas(:)  !! Mole fractions of each phase
      character(len=14), allocatable :: kinds_x(:)  !! Kinds of the main phases
      character(len=14) :: kind_w !! Kind of the reference phase
   end type MPEquilibriumState

contains

   type(MPEquilibriumState) function MPEquilibriumState_from_X(nc, np, z, kinds_x, kind_w, X)
      !! # `MPEquilibriumState_from_X`
      !!
      !! Function to create a `MPEquilibriumState` from the vector of variables
      !! `X` used in the multiphase generalized isopleths methods.
      integer, intent(in) :: nc !! Number of components
      integer, intent(in) :: np !! Number of phases - 1
      real(pr), intent(in) :: z(nc) !! Global composition
      character(len=14), intent(in) :: kinds_x(np) !! Kinds of the main phases
      character(len=14), intent(in) :: kind_w !! Kind of the reference phase
      real(pr), intent(in) :: X(nc*np+np+1+2) 
      !! Vector of variables.
      !! The first `nc*np` elements are the equilibrium constants for each
      !! phase and component, the next `np+1` elements are the mole fractions
      !! of the main phases (\(\beta^l\)), and the last two elements are 
      !! \(\lnP\) and \(\lnT\).

      integer :: l, lb, ub

      real(pr) :: K(np, nc)

      real(pr) :: x_l(np, nc), w(nc), betas(np+1), P, T

      do l=1,np
         lb = (l-1)*nc + 1
         ub = nc*np
         K(l, :) = exp(X(lb:ub))
      end do

      betas = X(np*nc+1 : np*nc+np+1)
      P = exp(X(nc*np+np+1 + 1))
      T = exp(X(nc*np+np+1 + 2))

      w = z/(matmul(betas(:np), K(:np, :)) + betas(np+1))
      do l=1,np
         x_l(l, :) = K(l, :) * w
      end do

      MPEquilibriumState_from_X = MPEquilibriumState(&
         z=z, &
         P=P, &
         T=T, &
         np=np, &
         x_l=x_l, &
         w=w, &
         betas=betas, &
         kinds_x=kinds_x, &
         kind_w=kind_w &
         )
   end function MPEquilibriumState_from_X

   subroutine write_EquilibriumState(eq, unit, iotype, v_list, iostat, iomsg)
      class(EquilibriumState), intent(in) :: eq
      integer, intent(in) :: unit
      character(*), intent(in) :: iotype
      integer, intent(in)  :: v_list(:)
      integer, intent(out) :: iostat
      character(*), intent(inout) :: iomsg

      character(*), parameter :: nl = new_line("G")

      write(unit, *) eq%kind, eq%T, eq%P, eq%beta, eq%x, eq%y

   end subroutine write_EquilibriumState
end module yaeos__equilibria_equilibrium_state
