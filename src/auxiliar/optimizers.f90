module yaeos__optimizers
   use yaeos__constants, only: pr
   implicit none

   type, abstract :: Optimizer
      logical :: verbose
      real(pr), allocatable :: parameter_step(:)
      real(pr) :: solver_tolerance = 1e-9_pr
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

module yaeos__optimizers_powell_wrap
   use yaeos__constants, only: pr
   use yaeos__optimizers, only: Optimizer, obj_func

   type, extends(Optimizer) :: PowellWrapper
      !! Wrapper derived type to optimize with the Powell method
   contains
      procedure :: optimize => powell_optimize
   end type PowellWrapper
contains
   subroutine powell_optimize(self, foo, X, F, data)
      use newuoa_module, only: newuoa
      class(PowellWrapper), intent(in out) :: self
      class(*), optional, intent(in out) :: data
      procedure(obj_func) :: foo
      real(pr), intent(in out) :: X(:)
      real(pr), intent(out) :: F

      real(pr) :: dx(size(x))

      integer :: n, npt 
      n = size(X)
      npt = (N+2 +(N+1)*(N+2)/2)/2

      if (allocated(self%parameter_step)) then
         dx = self%parameter_step
      else
         dx = X * 0.01_pr
      end if

      call newuoa(&
         n, npt, x, &
         maxval(abs(dx/10)), self%solver_tolerance, 0, int(1e9), wrap_foo &
      )
   contains
      subroutine wrap_foo(n, x, ff)
         integer  :: n
         real(pr) :: x(*)
         real(pr) :: ff
         real(pr) :: xx(n)
         xx = x(1:n)
         call foo(xx, FF, data=data)
         F = FF
      end subroutine wrap_foo
   end subroutine powell_optimize
end module