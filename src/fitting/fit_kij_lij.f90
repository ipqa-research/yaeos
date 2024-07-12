module yaeos__fitting_fit_kij_lij
   !! Binary interaction parameters fitting problem.
   use yaeos__fitting, only: pr, FittingProblem, ArModel
   implicit none

   integer, parameter :: nc=2
   type, extends(FittingProblem) :: FitKijLij
      !! # Binary Interaction Parameters of Cubic EoS fitting problem
      !! Fit the binary interaction parameters of a mixtures.
      !! 
      !! # Description
      !! Fitting setup for quadratic combining rules, it is possible to select
      !! which parameters will be optimized with the `fit_lij` and `fit_kij` 
      !! attributes.
      !! 
      !! # Examples
      !! 
      !! ## Fit the kij BIP
      !!
      !! ```fortran
      !!  type(CubicEoS) :: model ! Model to fit
      !!  type(FitKijLij) :: fitting_problem ! Fitting problem specification
      !!  type(EquilibriaState) :: exp_data(3)
      !!  real(pr) :: X(2) ! parameter variables
      !!  real(pr) :: error
      !!
      !!  ! <some procedure to define exp data>
      !! 
      !!  model = PengRobinson76(tc, pc, w) ! Model to fit
      !!
      !!  fitting_problem%exp_data = exp_data
      !!  fitting_problem%model = model
      !!  fitting_problem%fit_kij = .true.
      !!
      !!  X = 0 ! initial values == 0
      !!  err = optimize(X, fitting_problem)
      !! ```
      !!
      !! # References
      !!
      logical :: fit_lij = .false. !! Fit the \(l_{ij}\) parameter
      logical :: fit_kij = .false. !! Fit the \(k_{ij}\) parameter
   contains
      procedure :: get_model_from_X => model_from_X
   end type FitKijLij

contains

   subroutine model_from_X(problem, X)
      use yaeos, only: R, RKPR, PengRobinson78, ArModel, QMR, CubicEoS
      real(pr), intent(in) :: X(:)
      class(FitKijLij), intent(in out) :: problem

      real(pr) :: kij(nc, nc), lij(nc, nc)

      kij = 0
      kij(1, 2) = X(1)
      kij(2, 1) = kij(1, 2)

      lij = 0
      lij(1, 2) = X(2)
      lij(2, 1) = X(2)

      associate(model => problem%model)
      select type (model)
       class is (CubicEoS)
         associate (mr => model%mixrule)
            select type(mr)
             class is (QMR)
               if (problem%fit_kij) mr%k = kij
               if (problem%fit_lij) mr%l = lij
            end select
         end associate
      end select
      end associate
   end subroutine
end module yaeos__fitting_fit_kij_lij

