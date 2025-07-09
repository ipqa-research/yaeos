program main
   use yaeos
   implicit none
   type(CubicEoS) :: model
   type(QMRTD) :: mixrule

   integer, parameter :: nc=2
   real(pr) :: k0(nc, nc)=0, kinf(nc, nc)=0, Tref(nc, nc)=0

   real(pr) :: tc(nc), pc(nc), w(nc)

   real(pr) :: D, dDT, dDT2, dDi(nc), dDij(nc, nc), dDidT(nc)
   real(pr) :: ai(nc), daidt(nc), daidt2(nc)

   real(pr) :: n(nc), T, V, P, f1, f2, f3, f4, dx=1E-3
   integer :: i, j

   tc =  [126.2, 568.7]
   pc =  [33.98, 24.9]
   w  =  [9.01E-002,  0.492]

   k0(1, 2) = 0!-0.429315E+00
   k0(2, 1) = k0(1, 2)
   kinf(1, 2) = 0.485065E+01
   kinf(2, 1) = kinf(1, 2)

   Tref(1, 2) = Tc(1)
   Tref(2, 1) = Tref(1, 2)

   mixrule = QMRTD(k=kinf, k0=k0, Tref=Tref, l=0*k0)

   model = SoaveRedlichKwong(Tc, Pc, w, kij=kinf)
   n = [0.2, 0.8]
   T = 150
   

   print *, "==================================================================="
   print *, "NORMAL MR"
   call model%alpha%alpha(T/model%components%Tc, ai, daidt, daidt2)
   ai = ai*model%ac
   daidt = daidt*model%ac/model%components%Tc
   daidt2 = daidt2*model%ac/model%components%Tc**2
   call model%mixrule%Dmix(n, T, ai=ai, daidt=daidt, daidt2=daidt2, D=D, dDdT=dDT, dDdT2=dDT2, dDi=dDi, dDidT=dDidT, dDij=dDij)
   print *, D, dDT, dDT2
   print *, dDi
   print *, dDidT
   print *, dDij
   print *, "==================================================================="


   call model%set_mixrule(mixrule)

   call model%alpha%alpha(T/model%components%Tc, ai, daidt, daidt2)
   ai = ai*model%ac
   daidt = daidt*model%ac/model%components%Tc
   daidt2 = daidt2*model%ac/model%components%Tc**2


   call model%mixrule%Dmix(n, T, ai=ai, daidt=daidt, daidt2=daidt2, D=D, dDdT=dDT, dDdT2=dDT2, dDi=dDi, dDidT=dDidT, dDij=dDij)
   print *, D

   f2 = get_D(n, T+dx)
   f1 = get_D(n, T-dx)

   print *, (f2-f1)/(2*dx), dDT
   print *, (f2 - 2*D + f1)/dx**2, dDT2

   do i = 1, nc
      n = [0.2, 0.8]
      n(i) = n(i) + dx
      f2 = get_D(n, T)
      n(i) = n(i) - 2*dx
      f1 = get_D(n, T)
      n(i) = n(i) + dx
      print *, (f2-f1)/(2*dx), dDi(i)
   end do

   dx = 0.01

   f1 = get_D(n + [dx,   dx], T)
   f2 = get_D(n + [dx,  -dx], T)
   f3 = get_D(n + [-dx,  dx], T)
   f4 = get_D(n + [-dx, -dx], T)
   print *, (f1 - f2 - f3 + f4)/(4*dx**2), dDij(1, 2), dDij(2, 1)

   dx = 0.01
   f1 = get_D(n + [dx,   0._pr], T)
   f2 = get_D(n - [dx,   0._pr], T)
   print *, (f1 - 2*D + f2)/(dx**2), dDij(1, 1)
   
   f1 = get_D(n + [0._pr,   dx], T)
   f2 = get_D(n - [0._pr,   dx], T)
   print *, (f1 - 2*D + f2)/(dx**2), dDij(2, 2)


contains
   real(pr) function get_D(n, T) result(D)
      real(pr) :: n(nc), T

      real(pr) :: ai(nc), daidt(nc), daidt2(nc)
      real(pr) :: dDT, dDT2, dDi(nc), dDij(nc, nc), dDidT(nc)

      call model%alpha%alpha(T/model%components%Tc, ai, daidt, daidt2)
      ai = ai*model%ac
      daidt = daidt*model%ac/model%components%Tc
      daidt2 = daidt2*model%ac/model%components%Tc**2
      call model%mixrule%Dmix(n, T, ai=ai, daidt=daidt, daidt2=daidt2, D=D, dDdT=dDT, dDdT2=dDT2, dDi=dDi, dDidT=dDidT, dDij=dDij)
   end function get_D
end program main
