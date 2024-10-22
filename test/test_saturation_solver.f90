program main
   use yaeos
   use yaeos__math, only: solve_system
   use yaeos__m_s_sp, only: saturation_F, saturation_TP, solve_TP, solve_VxVyT
   use fixtures_models, only: binary_PR76
   use fortime, only: timer

   implicit none

   class(ArModel), allocatable :: model
   type(EquilibriumState) :: sat

   integer, parameter :: nc = 2, nf = nc + 3

   real(pr) :: n(nc), T, P, Vy, Vz, X(nf), S, F(nf), df(nf,nf)
   real(pr) :: dx(nf), dFnum(nf, nf), Fdx(nf), dftmp(nf,nf)
   real(pr) :: Px, Py
   integer :: i, ns
   type(timer) :: tim

   model = binary_PR76()

   n = [0.5, 0.5]

   T = 200._pr

   sat = saturation_pressure(model, n, T=T, kind="bubble")

   print *, "saturation_pressure"
   print *, sat%iters, sat%x, sat%y, sat%P

   print *, "solve yVxVy"

   ! call numdiff

   call VxVy
   call TP

contains

   subroutine numdiff
      Vz = sat%Vx! *0.9
      Vy = sat%Vy
      call model%pressure(n=n, V=Vz, T=T, P=P)

      X(:nc) = log(sat%y/sat%x)
      X(nc+1) = Vz
      X(nc+2) = Vy
      X(nc+3) = T
      S = T
      ns = nc+3
      dx = 0
      call saturation_F(model, n, X, ns, S, F, dF)

      print *, F
      print *, "numdiff"
      do i=1,nf
         dx = 0
         dx(i) = 1e-9

         call saturation_F(model, n, X+dx, ns, S, Fdx, dFtmp)
         dFnum(:, i) = (Fdx - F) / dx(i)

         print *, dfnum(:, i)
         print *, df(:, i)
         print *, "=========================="
      end do
   end subroutine numdiff

   subroutine VxVy
      real(pr) :: X(nc+3), S, tol=1e-10
      integer :: ns, its
      X(1) = 1.0
      X(2) = 0.5
      X(3) = log(0.1)
      X(4) = log(0.7)
      X(5) = T

      ns = 5
      S = X(ns)

      call solve_VxVyT(model, n, X, ns, S, tol, 100, its)
      print *, its, n*exp(X(:nc)), exp(X(nc+1)), exp(X(nc+2)), X(5)
   end subroutine VxVy

   subroutine TP
      real(pr) :: X(nc+2), S, tol=1e-10
      integer :: its, ns
      print *, "solve_TP"
      X(1) = log(1.0)
      X(2) = log(0.1)
      X(3) = log(T)
      X(4) = log(200.0)
      ns = 3

      S = X(ns)

      call solve_TP(model, "bubble", n, X, ns, S, tol, 100, its)
      print *, its, n*exp(X(:nc)), exp(X(nc+1:))
   end subroutine TP
end program main
