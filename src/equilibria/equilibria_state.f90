module yaeos__equilibria_equilibrium_state
   !! # Equilibrium state types
   !!
   !! # Description
   !! Derived types used to store the results of phase-equilibrium
   !! calculations. [[EquilibriumState]] holds a two-phase state and
   !! [[MPEquilibriumState]] holds a multiphase state.
   !! The values use the library's default units: pressure in bar, temperature
   !! in kelvin, volume in litres per mole, and compositions as mole fractions.
   !!
   !! # Examples
   !! ```fortran
   !! result = flash(model, z, T, p_spec=P)
   !! print *, result%kind, result%beta, result%x, result%y
   !! ```
   use yaeos__constants, only: pr
   implicit none

   type :: EquilibriumState
      !! # EquilibriumState
      !! Result of a two-phase flash or saturation-point calculation.
      !!
      !! # Description
      !! Stores the state kind, phase compositions, volumes, conditions, phase
      !! fraction, and iteration count. `x` denotes the heavy phase and `y` the
      !! light phase; `beta` is the light-phase mole fraction.
      !!
      !! # Examples
      !! ```fortran
      !! print *, result%P, result%T, result%beta
      !! ```
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
      !! # MPEquilibriumState
      !! Result of a multiphase equilibrium calculation.
      !!
      !! # Description
      !! Holds the feed, conditions, main-phase compositions and kinds, plus
      !! a reference phase used in the generalized multiphase formulation.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! integer :: phase
      !! result = pt_mp_flash(model, z, P, T)
      !! do phase = 1, result%np
      !!    print *, result%kinds_x(phase), result%x_l(phase, :)
      !! end do
      !! ```
      real(pr), allocatable :: z(:) !! Global composition
      real(pr) :: P !! Pressure [bar]
      real(pr) :: T !! Temperature [K]
      integer :: np !! Number of main phases (excluding the reference phase)
      real(pr), allocatable :: x_l(:,:) !! Mole fractions of the main phases, shape `(np, nc)`
      real(pr), allocatable :: w(:) !! Mole fractions of the reference phase
      real(pr), allocatable :: betas(:) !! Phase fractions: `np` main phases followed by the reference phase
      character(len=14), allocatable :: kinds_x(:) !! Kinds of the main phases
      character(len=14) :: kind_w !! Kind of the reference phase
   end type MPEquilibriumState

contains

   type(MPEquilibriumState) function MPEquilibriumState_from_X(nc, np, z, kinds_x, kind_w, X)
      !! # `MPEquilibriumState_from_X`
      !!
      !! Function to create a `MPEquilibriumState` from the vector of variables
      !! `X` used in the multiphase generalized isopleths methods.
      integer, intent(in) :: nc !! Number of components
      integer, intent(in) :: np !! Number of main phases (total phases - 1)
      real(pr), intent(in) :: z(nc) !! Global composition
      character(len=14), intent(in) :: kinds_x(np) !! Kinds of the main phases
      character(len=14), intent(in) :: kind_w !! Kind of the reference phase
      real(pr), intent(in) :: X(nc*np+np+1+2)
      !! Vector of variables.
      !! The first `nc*np` elements are the equilibrium constants for each
      !! phase and component, the next `np+1` elements are the mole fractions
      !! of the main phases (\(\beta^l\)), and the last two elements are 
      !! \(\lnP\) and \(\lnT\).

      integer :: l !! Loop/phase index.
      integer :: lb
      integer :: ub

      real(pr) :: K(np, nc)

      real(pr) :: x_l(np, nc) !! Main-phase compositions.
      real(pr) :: w(nc) !! Reference/trial phase composition.
      real(pr) :: betas(np+1) !! Phase fractions of the main phases followed by the reference phase.
      real(pr) :: P
      real(pr) :: T

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
      !! # write_EquilibriumState
      !! Format a two-phase equilibrium state for formatted derived-type I/O.
      !!
      !! # Description
      !! Writes the phase kind, temperature, pressure, phase fraction, and
      !! phase compositions to the requested Fortran I/O unit.
      class(EquilibriumState), intent(in) :: eq !! State to serialize.
      integer, intent(in) :: unit !! Output unit number.
      character(*), intent(in) :: iotype !! Fortran formatted-I/O mode string.
      integer, intent(in)  :: v_list(:) !! Format control list supplied by the runtime.
      integer, intent(out) :: iostat !! I/O status returned to the runtime.
      character(*), intent(inout) :: iomsg !! I/O error message buffer.

      character(*), parameter :: nl = new_line("G")

      write(unit, *) eq%kind, eq%T, eq%P, eq%beta, eq%x, eq%y

   end subroutine write_EquilibriumState
end module yaeos__equilibria_equilibrium_state
