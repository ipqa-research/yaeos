module yaeos__fitting
   use yaeos__constants, only: pr
   use yaeos__models, only: ArModel
   use yaeos__equilibria, only: &
      EquilibriumState, saturation_pressure, saturation_temperature, flash
   use yaeos__optimizers, only: Optimizer, obj_func
   use ieee_arithmetic, only: isnan => ieee_is_nan
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
      !! Residual Helmholtz Model to fit
      type(EquilibriumState), allocatable :: experimental_points(:)
      !! Experimental points to fit
      logical :: verbose = .false.
      !! If true log the fitting process
   contains
      procedure(model_from_X), deferred :: get_model_from_X
   end type FittingProblem

   abstract interface
      subroutine model_from_X(problem, X)
         !! Function that returns a setted model from the parameters vector
         import ArModel, FittingProblem, pr
         class(FittingProblem), intent(in out) :: problem
         !! Fitting problem to optimize
         real(pr), intent(in) :: X(:)
         !! Vector of parameters to fit
      end subroutine model_from_X
   end interface

contains

   real(pr) function optimize(X, opt, data) result(y)
      real(pr), intent(in out) :: X(:)
      !! Vector of parameters to fit
      class(Optimizer), intent(in out) :: opt
      !! Optimizer object, bsaed on the `Optimizer` class from
      !! `yaeos__optimizers`
      class(FittingProblem), optional, intent(in out) :: data
      !! Fitting problem to optimize
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
      use yaeos__equilibria_saturation_points, only: max_iterations
      real(pr), intent(in) :: X(:) !! Vector of parameters
      real(pr), intent(out) :: Fobj !! Objective function
      real(pr), optional, intent(out) :: dF(:)
      !! Gradient of the objective function, only exists to be consistent
      !! with the `Optimizer` class API
      class(*), optional, intent(in out):: func_data

      type(EquilibriumState) :: model_point !! Each solved point
      type(EquilibriumState) :: exp_point

      integer :: i
      logical :: pt_converged
      integer :: n_nconv

      if (present(dF)) error stop 1

      select type(func_data)
       class is (FittingProblem)
         block
            real(pr) :: fobjs(size(func_data%experimental_points))

            ! Update the problem model to the new vector of parameters
            call func_data%get_model_from_X(X)

            n_nconv = 0
            fobj = 0
            ! Calculate each point and  calculate its error.
            ! if at some point there is a NaN value, assign a big number

            do i=1, size(func_data%experimental_points)
               associate( model => func_data%model )

                  exp_point = func_data%experimental_points(i)

                  select case(exp_point%kind)
                   case("bubble")
                     model_point = saturation_pressure(&
                        model, exp_point%x, exp_point%t, kind="bubble", &
                        p0=exp_point%p, y0=exp_point%y)
                   case("dew")
                     model_point = saturation_pressure(&
                        model, exp_point%y, exp_point%t, kind="dew", &
                        p0=exp_point%p, y0=exp_point%x)
                   case("liquid-liquid")
                     model_point = saturation_pressure(&
                        model, exp_point%x, exp_point%t, kind="liquid-liquid", &
                        p0=exp_point%p, y0=exp_point%y)
                  end select


                  if (model_point%iters > max_iterations .or. isnan(model_point%P)) then
                     ! If the point did not converge, calculate the phase
                     ! envelope and get the closest point
                     call calc_pt_envel(model, exp_point, model_point, pt_converged)

                     fobjs(i) = sq_error(exp_point%p, model_point%p) + sq_error(exp_point%T, model_point%T)

                     if (.not. pt_converged) then
                        n_nconv = n_nconv + 1
                        fobjs(i) = maxval(func_data%experimental_points%P)
                     end if

                  else
                     ! Calculate the error
                     fobjs(i) = sq_error(exp_point%p, model_point%p)
                  end if
                  ! print *, exp_point%T, model_point%T, exp_point%P, model_point%P, fobjs(i)

                  if (isnan(fobjs(i))) fobjs(i) = 1e25
               end associate
            end do


            fobjs = fobjs * func_data%experimental_points%P**2
            fobj = sum(fobjs)/size(fobjs)

            if (func_data%verbose) then
               write(*, "(I3,2x(E20.9),2x, '[',*(E15.6,','))", advance="no") n_nconv, fobj, X
               write(*, "(']', /, '==========================================')")
            end if

         end block
      end select

   end subroutine error_function

   subroutine calc_pt_envel(model, exp_point, model_point, converged)
      use yaeos__math, only: sq_error, interpol
      use yaeos__equilibria, only: PTEnvel2, EquilibriumState, pt_envelope_2ph, k_wilson
      use yaeos__equilibria_saturation_points, only: max_iterations
      class(ArModel), intent(in) :: model
      type(EquilibriumState), intent(in) :: exp_point
      type(EquilibriumState), intent(out) :: model_point
      logical, intent(out) :: converged

      type(PTEnvel2) :: env
      type(EquilibriumState) :: init
      integer :: pos
      integer, allocatable :: msk(:)

      real(pr) :: z(size(exp_point%x))

      select case(exp_point%kind)
       case("bubble")
         z = exp_point%x
       case("dew")
         z = exp_point%y
       case("liquid-liquid")
         z = exp_point%x
      end select

      init = saturation_temperature(model, z, P=1.0_pr, kind="dew", T0=500._pr)
      if (init%iters > max_iterations) then
         init = saturation_pressure(model, z, T=150._pr, kind="bubble")
      end if
      env = pt_envelope_2ph(model, z, init, points=2000)

      if (size(env%points) > 50) then
         converged = .true.
         pos = minloc(abs(env%points%T - exp_point%T), dim=1)

         ! Find closest point to the experimental point
         pos = minloc((env%points%T - exp_point%T)**2 + (env%points%P - exp_point%P)**2, dim=1)
         model_point = env%points(pos)

         ! Interpolation
         !model_point%P = interpol(env%points(pos)%T, env%points(pos+1)%T, &
         !                         env%points(pos+1)%P, env%points(pos+1)%P, exp_point%T)

         ! if (abs(model_point%T - exp_point%T) > 5) then
         !    converged = .false.
         ! end if
      else
         converged = .false.
      end if
   end subroutine calc_pt_envel
end module yaeos__fitting
