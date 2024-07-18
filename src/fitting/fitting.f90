module yaeos__fitting
   use yaeos__constants, only: pr
   use yaeos__models, only: ArModel, CubicEoS
   use yaeos__equilibria, only: &
      EquilibriumState, saturation_pressure, saturation_temperature, flash
   use forbear, only: bar_object
   implicit none

   type, abstract :: FittingProblem
      !! # Fitting problem setting
      !!
      !! # Description
      !! This derived type holds all the relevant information for a parameter
      !! optimization problem. It keeps the base model structure that will be
      !! optimized and a procedure `get_model_from_X` that should reconstruct
      !! the model with the desired parameters to optimize.
      real(pr) :: solver_tolerance = 1e-9_pr
      real(pr), allocatable :: parameter_step(:)

      class(ArModel), allocatable :: model

      type(EquilibriumState), allocatable :: experimental_points(:)
      logical :: verbose = .false.
   contains
      procedure(model_from_X), deferred :: get_model_from_X
   end type FittingProblem

   abstract interface
      subroutine model_from_X(problem, X)
         !! Function that returns a setted model from the parameters vector
         import ArModel, FittingProblem, pr
         class(FittingProblem), intent(in out) :: problem
         real(pr), intent(in) :: X(:)
      end subroutine model_from_X
   end interface

   type(bar_object), private :: bar
   integer, private :: count

contains

   real(pr) function optimize(X, func_data) result(y)
      use nlopt_wrap, only: create, destroy, nlopt_opt, nlopt_algorithm_enum
      use nlopt_callback, only: nlopt_func, create_nlopt_func

      real(pr), intent(in out) :: X(:) !! Vector of parameters to fit
      class(FittingProblem) :: func_data !! Parametrization details

      real(pr) :: dx(size(X))

      type(nlopt_opt) :: opt !! Optimizer
      type(nlopt_func) :: f !! Function to optimize
      integer :: stat
      integer :: i

      count = 0
      call bar%initialize(&
         prefix_string='Fitting... ',&
         width=1, spinner_string='⠋', spinner_color_fg='blue', &
         min_value=0._pr, max_value=100._pr &
         )
      call bar%start

      ! opt = nlopt_opt(nlopt_algorithm_enum%LN_NELDERMEAD, size(X))
      ! opt = nlopt_opt(nlopt_algorithm_enum%LN_BOBYQA, size(X))
      ! opt = nlopt_opt(nlopt_algorithm_enum%LN_NEWUOA, size(X))
      opt = nlopt_opt(nlopt_algorithm_enum%LN_PRAXIS, size(X))

      f = create_nlopt_func(fobj, f_data=func_data)

      if (.not. allocated(func_data%parameter_step)) then
         dx = [(0.05_pr, i=1, size(X))]
      else
         dx = func_data%parameter_step
      end if

      call opt%set_ftol_rel(func_data%solver_tolerance)

      call opt%set_initial_step(dx)
      call opt%set_min_objective(f)
      call opt%set_xtol_abs(dx/100)
      call opt%optimize(x, y, stat)
      call bar%destroy
      call destroy(opt)
   end function optimize

   real(pr) function fobj(x, gradient, func_data)
      !! # Objective function to fit phase-equilibria points.
      !!
      !! # Description
      !! ...
      !!
      !! # Examples
      !!
      !! ```fortran
      !!  ...
      !! ```
      use yaeos__math, only: sq_error
      real(pr), intent(in) :: x(:)
      real(pr), optional, intent(in out) :: gradient(:)
      class(*), optional, intent(in) :: func_data

      integer :: i

      real(pr) :: p_exp, t_exp

      select type(func_data)
       class is(FittingProblem)
         fobj = error_function(X, func_data)
         if (func_data%verbose) then
            call bar%update(current=real(count,pr)/(count + 100))
            write(*, "(E15.4, 2x)", advance="no") fobj
         end if
      end select
      write(2, "(*(E15.4,2x))") X, fobj
      write(1, "(/)")

      count = count + 1
   end function fobj

   real(pr) function error_function(X, func_data) result(fobj)
      !! # `error_function`
      !! Error function for phase-equilibria optimization. Using two-phase
      !! points and an error function of:
      !! 
      !! \[
      !!   FO =   \sum_i (\frac{P_i^{exp} - P_i^{calc}}{P_i^{exp}})^2
      !!        + \sum_i (y_i^{exp} - y_i^{calc})**2
      !!        + \sum_i (x_i^{exp} - x_i^{calc})**2
      !! \]
      use yaeos__math, only: sq_error
      real(pr), intent(in) :: X(:)
      class(FittingProblem) :: func_data

      type(EquilibriumState) :: model_point !! Each solved point
      type(EquilibriumState) :: exp_point

      integer :: i

      ! Update the problem model to the new vector of parameters
      call func_data%get_model_from_X(X)

      fobj = 0

      associate( model => func_data%model )
         
         ! Calculate each point and  calculate its error.
         ! if at some point there is a NaN value, assign a big number and
         ! exit
         do i=1, size(func_data%experimental_points)
            exp_point = func_data%experimental_points(i)

            select case(exp_point%kind)
             case("bubble")
               model_point = saturation_pressure(&
                  model, exp_point%x, exp_point%t, kind="bubble", &
                  p0=exp_point%p, y0=exp_point%y &
                  )
             case("dew")
               model_point = saturation_pressure(&
                  model, exp_point%y, exp_point%t, kind="dew", &
                  p0=exp_point%p, y0=exp_point%x &
                  )
             case("liquid-liquid")
               model_point = saturation_pressure(&
                  model, exp_point%x, exp_point%t, kind="liquid-liquid", &
                  p0=exp_point%p, y0=exp_point%y &
                  )
            end select

            fobj = fobj + sq_error(exp_point%p, model_point%p)
            fobj = fobj + maxval(sq_error(exp_point%y, model_point%y))
            fobj = fobj + maxval(sq_error(exp_point%x, model_point%x))
            write(1, *) exp_point, model_point

            if(isnan(fobj)) then
               ! fobj = 1e6
               exit
            end if
         end do
      end associate
   end function error_function
end module yaeos__fitting
