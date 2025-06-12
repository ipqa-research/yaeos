program main
   use yaeos
   use yaeos__math, only: solve_system
   use yaeos__m_s_sp, only: saturation_F, saturation_TP, solve_TP, solve_VxVyT
   use fixtures_models, only: binary_PR76
   use fortime, only: timer

   implicit none

   class(ArModel), allocatable :: model
   type(EquilibriumState) :: sat
   type(PTEnvel2) :: env

   integer, parameter :: nc = 2, nf = nc + 4

   real(pr) :: n(nc), T, P, Vy, Vz, dPdVy, dPdVz, X(nf), S, F(nf), df(nf,nf)
   real(pr) :: dx(nf), dFnum(nf, nf), Fdx(nf), dftmp(nf,nf)
   real(pr) :: Px, Py
   integer :: i, ns
   type(timer) :: tim
   real(pr) :: et,st, t1, t2

   model = binary_PR76()

   n = [0.5, 0.5]

   T = 200._pr

   sat = saturation_pressure(model, n, T=T, kind="bubble")
   sat = saturation_temperature(model, n, P=20._pr, kind="bubble")


   print *, "saturation_pressure"
   print *, sat%iters, sat%x, sat%y, sat%T, sat%P
   
   ! call numdiff

   call cpu_time(st)
   print *, "solve yVxVy"
   call VxVy
   call cpu_time(et)
   t1 = et-st
   print *, (et-st)*1e6, "us"
   
   call cpu_time(st)
   call TP
   call cpu_time(et)
   t2 = et-st
   print *, (et-st)*1e6, "us"

   print *, t2/t1

contains

   subroutine numdiff
      Vz = sat%Vx! *0.9
      Vy = sat%Vy
      call model%pressure(n=n, V=Vz, T=T, P=P)

      X(:nc) = log(sat%y/sat%x)
      X(nc+1) = log(Vz)
      X(nc+2) = log(Vy)
      X(nc+3) = log(T)
      X(nc+4) = log(P)
      S = T
      ns = nc+3
      dx = 0
      call saturation_F(model, n, X, ns, S, F, dF, dPdVz, dPdVy)

      print *, F
      print *, "numdiff"
      do i=1,nf
         dx = 0
         dx(i) = 1e-9

         call saturation_F(model, n, X+dx, ns, S, Fdx, dFtmp, dPdVz, dPdVy)
         dFnum(:, i) = (Fdx - F) / dx(i)

         print *, dfnum(:, i)
         print *, df(:, i)
         print *, "=========================="
      end do
   end subroutine numdiff

   subroutine VxVy
      real(pr) :: X(nc+4), S, tol=1e-10
      integer :: ns, its
      X(1) = 1.0
      X(2) = 0.5
      X(3) = log(0.1)
      X(4) = log(0.7)
      X(5) = log(T)
      X(6) = log(20.0)

      ns = 6
      S = X(ns)

      print *, "solveVxVyPT"

      call solve_VxVyT(model, n, X, ns, S, tol, 1000, its)
      print *, its, n*exp(X(:nc)), exp(X(5)), exp(X(6))
   end subroutine VxVy

   subroutine TP
      real(pr) :: X(nc+2), S, tol=1e-10

      real(pr) :: lnphi_z(nc), lnphi_y(nc)

      real(pr) :: Vz, Vy, y(nc), Pz, Py
      integer :: its, ns
      print *, "solve_TP"
      X(1) = log(1.0)
      X(2) = log(0.1)
      X(3) = log(230.)
      X(4) = log(20.0)
      ns = 4

      S = X(ns)

      call solve_TP(model, "bubble", n, X, ns, S, tol, 1000, its)
      y = n * exp(X(:nc))
      print *, its, y, exp(X(nc+1:))
   end subroutine TP
end program main
