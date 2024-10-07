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
      !!  type(EquilibriumState) :: exp_data(3)
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
      logical :: fit_kij_exp_t = .false. !! Fit \(k_{ij}(T)\)
   contains
      procedure :: get_model_from_X => model_from_X
   end type FitKijLij

contains

   subroutine model_from_X(problem, X)
      use yaeos, only: pr, CubicEoS, QMR, QMR_RKPR
      use yaeos__models_ar_cubic_quadratic_mixing, only: QMRTD
      real(pr), intent(in) :: X(:)
      class(FitKijLij), intent(in out) :: problem

      real(pr) :: kij(nc, nc), lij(nc, nc), Tref(nc,nc)

      if (size(X) > 3) error stop 1

      kij = 0
      lij = 0

      if (problem%fit_kij .and. problem%fit_lij) then
         kij(1, 2) = X(1)
         kij(2, 1) = kij(1, 2)
         lij(1, 2) = X(2)
         lij(2, 1) = lij(1, 2)
      else if (problem%fit_kij) then
         kij = 0
         kij(1, 2) = X(1)
         kij(2, 1) = kij(1, 2)
      else if (problem%fit_lij) then
         lij = 0
         lij(1, 2) = X(1)
         lij(2, 1) = lij(2, 1)
      end if

      associate(model => problem%model)
         select type (model)
          class is (CubicEoS)
            associate (mr => model%mixrule)
               select type(mr)
                type is (QMR)
                  if (problem%fit_kij) mr%k = kij
                  if (problem%fit_lij) mr%l = lij
                type is (QMR_RKPR)
                  if (problem%fit_kij) mr%k = kij
                  if (problem%fit_lij) mr%l = lij
                  ! call model%set_delta1(X(2:))
                type is (QMRTD)
                  if (problem%fit_kij) mr%k = kij   ! kinf
                  if (problem%fit_kij) mr%k0  = lij ! k0
               end select

            end associate
         end select
      end associate
   end subroutine model_from_X
end module yaeos__fitting_fit_kij_lij

