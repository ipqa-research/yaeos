module yaeos__fitting
   use yaeos__constants, only: pr
   use yaeos__models, only: ArModel, CubicEoS
   use yaeos__equilibria, only: &
      EquilibriumState, saturation_pressure, saturation_temperature, flash
   use yaeos__optimizers, only: Optimizer, obj_func
   implicit none

   type, abstract :: FittingProblem
      !! # Fitting problem setting
      !!
      !! # Description
      !! This derived type holds all the relevant information for a parameter
      !! optimization problem. It keeps the base model structure that will be
      !! optimized and a procedure `get_model_from_X` that should reconstruct
      !! the model with the desired parameters to optimize.
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

contains

   real(pr) function optimize(X, opt, data) result(y)
      real(pr), intent(in out) :: X(:) !! Vector of parameters to fit
      class(Optimizer), intent(in out) :: opt
      class(FittingProblem), intent(in out) :: data

      call opt%optimize(error_function, X, y, data)
   end function optimize

   subroutine error_function(X, Fobj, dF, func_data)
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
      real(pr), intent(out) :: Fobj
      real(pr), optional, intent(out) :: dF(:)
      class(*), intent(in out):: func_data

      type(EquilibriumState) :: model_point !! Each solved point
      type(EquilibriumState) :: exp_point

      integer :: i

      if (present(dF)) error stop 1

      select type(func_data)
       class is (FittingProblem)
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
                  fobj = 1e6
                  exit
               end if
            end do
         end associate
      end select
   end subroutine error_function
end module yaeos__fitting
