!! Example of adding a new optimizer to the library. In this case is the
!! wrapper for the NLOPT library. NLopt is a free/open-source library for
!! nonlinear optimization, providing a common interface for a number of
!! different optimization algorithms available in the library.
module yaoes__optimizers_nlopt_wrap
   use yaeos__constants, only: pr
   use yaeos__optimizers, only: Optimizer, obj_func
   
   type, extends(Optimizer) :: NLOPTWrapper
      !! Wrapper derived type to optimize with the NLOP library
   contains
      procedure :: optimize => nlopt_optimize
   end type NLOPTWrapper
contains
   subroutine nlopt_optimize(self, foo, X, F, data)
      use nlopt_wrap, only: create, destroy, nlopt_opt, nlopt_algorithm_enum
      use nlopt_callback, only: nlopt_func, create_nlopt_func
      class(NLOPTWrapper), intent(in out) :: self
      class(*), optional, intent(in out) :: data
      procedure(obj_func) :: foo
      real(pr), intent(in out) :: X(:)
      real(pr), intent(out) :: F

      real(pr) :: dx(size(x))
      integer :: i, stat

      type(nlopt_opt) :: opt
      type(nlopt_func) :: fun

      fun = create_nlopt_func(nlopt_func_wrapper, f_data=data)
      opt = nlopt_opt(nlopt_algorithm_enum%LN_PRAXIS, size(X))

      call opt%set_ftol_rel(self%solver_tolerance)
      
      if (allocated(self%parameter_step)) then
         dx = self%parameter_step
      else
         dx = 0.01_pr
      end if

      call opt%set_min_objective(fun)
      call opt%set_xtol_abs(dx/100)
      call opt%optimize(x, F, stat)
      call destroy(opt)
   contains
      real(pr) function nlopt_func_wrapper(xx, gradient, func_data) result(fobj)
         real(pr), intent(in) :: xX(:)
         real(pr), optional, intent(in out) :: gradient(:)
         class(*), optional, intent(in) :: func_data
         call foo(xx, fobj, gradient, data)
      end function nlopt_func_wrapper
   end subroutine nlopt_optimize
end module

