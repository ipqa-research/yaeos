module yaeos__math_nonlinearsolvers
   use yaeos__constants, only: pr
   use yaeos__math_linalg, only: solve_system
   implicit none

   interface to_solve
      subroutine to_solve(X, F, J)
         import pr
         real(pr), intent(in) :: X(:)
         real(pr), intent(out) :: F(:)
         real(pr), intent(out) :: J(:, :)
      end subroutine to_solve
   end interface

contains

   subroutine newton(sub, x, tol, max_its, its)
      procedure(to_solve) :: sub
      real(pr), intent(in out) :: x(:)
      real(pr), intent(in) :: tol
      integer, intent(in) :: max_its
      integer, intent(out) :: its

      real(pr) :: F(size(X)), J(size(X), size(X)), dX(size(X))

      do its = 1, max_its
         call sub(x, F, J)
         if (maxval(abs(F)) < tol) exit
         dX = solve_system(J, -F)
         X = X + dX
         write(1, *) its, maxval(abs(F))
      end do
   end subroutine newton

   subroutine homotopy(sub, x, tol, max_its, its)
      procedure(to_solve) :: sub
      real(pr), intent(in out) :: x(:)
      real(pr), intent(in) :: tol
      integer, intent(in) :: max_its
      integer, intent(out) :: its

      real(pr) :: dX(size(x))
      real(pr) :: F(size(X)), dF(size(X), size(X))
      real(pr) :: G(size(X)), dG(size(X), size(X))
      real(pr) :: H(size(X)), dH(size(X), size(X))
      real(pr) :: t, dt
      real(pr) :: X0(size(X))
      integer :: i, nt

      nt = 5
      dt = 1.0_pr/real(nt, pr)


      do i=0,nt
         t = dt*i
         X0 = X
         call newton(wrap, X, tol, 15, its)
         call wrap(X, H, dH)
         print *,  t, maxval(abs(H))
         print *, " ====== "
      end do
   contains
      subroutine wrap(X, H, dH)
         real(pr), intent(in) :: X(:)
         real(pr), intent(out) :: H(:)
         real(pr), intent(out) :: dH(:,:)

         real(pr) :: F(size(X)), dF(size(X), size(X))
         real(pr) :: G(size(X)), dG(size(X), size(X))
         real(pr) :: diag(size(X), size(X)), id(size(X), size(X))
         integer :: i

         G = 0
         G = X - X0
         do i=1,size(X)
            dG(i, i) = 1
         end do

         call sub(X, F, dF)

         print *, G

         H = t * F + (1-t) * G
         dH = t * dF + (1-t) * dG

         ! H = (t * F + (1-t) * G)**2

         ! do i=1,size(X)
         !    diag(i, i) = H(i)
         ! end do
         ! dH = 2 * diag * (t * dF + (1-t) * dG)
      end subroutine wrap
   end subroutine homotopy
end module yaeos__math_nonlinearsolvers
