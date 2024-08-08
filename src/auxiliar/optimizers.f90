module yaeos__optimizers
   use yaeos__constants, only: pr
   implicit none

   type, abstract :: Optimizer
   contains
      procedure(abs_optimize), deferred :: optimize
   end type

   abstract interface
      subroutine obj_func(X, F, dF)
         import pr
         real(pr), intent(in) :: X(:)
         real(pr), intent(out) :: F
         real(pr), optional, intent(out) :: dF(:)
      end subroutine
   end interface

   abstract interface
      subroutine abs_optimize(self, foo, X, F)
         import pr, obj_func, Optimizer
         class(Optimizer), intent(in out) :: self
         procedure(obj_func) :: foo
         real(pr), intent(in) :: X(:)
         real(pr), intent(out) :: F
      end subroutine
   end interface
end module