module yaeos__optimizers
   use yaeos__constants, only: pr
   implicit none

   type, abstract :: Optimizer
   contains
      procedure(abs_optimize), deferred :: optimize
   end type

   abstract interface
      subroutine obj_func(X, F, dF, data)
         import pr
         real(pr), intent(in) :: X(:)
         real(pr), intent(out) :: F
         real(pr), optional, intent(out) :: dF(:)
         class(*), optional, intent(in out) :: data
      end subroutine
   end interface

   abstract interface
      subroutine abs_optimize(self, foo, X, F, data)
         import pr, obj_func, Optimizer
         class(Optimizer), intent(in out) :: self
         procedure(obj_func) :: foo
         real(pr), intent(in out) :: X(:)
         real(pr), intent(out) :: F
         class(*), optional, intent(in out) :: data
      end subroutine
   end interface
end module


module yaoes__optimizers_nlopt_wrap
   use yaeos__constants, only: pr
   use yaeos__optimizers
   
   type, extends(Optimizer) :: NLOPTWrapper
      !! Wrapper derived type to optimize with the NLOP library
      real(pr) :: solver_tolerance = 1e-9_pr
      real(pr), allocatable :: parameter_step(:)
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