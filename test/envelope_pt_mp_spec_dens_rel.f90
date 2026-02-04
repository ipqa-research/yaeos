program main
   use testing_aux, only: assert, test_title
   use yaeos
   use fixtures_models, only: binary_PR76
   implicit none
   type(CubicEoS) :: model
   type(EquilibriumState) :: sat

   integer, parameter :: nc=2, np=2
   integer, parameter :: psize = np*nc + np + 2
   real(pr) :: z(nc), x_l0(np, 2), w0(nc), p0, t0, betas0(np)
   real(pr) :: X(nc*np + np + 2), Vnv, Vw
   integer :: ns
   real(pr) :: S
   real(pr) :: maxerr

   print *, test_title("Specify ns < 0 in PT multi-phase envelope")
   model = binary_PR76()

   t0 = 150

   z = [0.5, 0.5]

   sat = saturation_pressure(model, n=z, T=t0, kind="bubble", p0=1._pr)

   Vnv = sat%Vx
   Vw = sat%Vy
   S = log(Vnv/Vw)
   ns = -1

   X(1:nc) = log(z/sat%y)
   X(nc+1:nc*np) = log(sat%y/z)
   X(np*nc + 1) = 0.3_pr
   X(np*nc + 2) = 0.7_pr
   X(np*nc + np + 1) = log(sat%P)
   X(np*nc + np + 2) = log(sat%T)

   call numdiff(maxerr)
   call assert(maxerr < 1e-4_pr, "Max error in numerical differentiation is too large")

contains
   subroutine numdiff(maxerr)
      use yaeos__equilibria_boundaries_phase_envelopes_mp, only: pt_F_NP
      real(pr), intent(out) :: maxerr
      real(pr) :: dfnum(psize, psize), F1(psize), F2(psize), F(psize)
      real(pr) :: tmp(psize, psize)
      real(pr) :: XdX(psize)
      real(pr) :: eps = 1e-3, dx
      real(pr) :: Vl(np), Vw
      character(len=*), parameter :: fmt="(*(E10.2,2x))"
      character(len=14) :: kinds_x(np)
      character(len=14) :: kind_w
      real(pr) :: df(psize, psize)


      integer :: i, j, loc(2)

      kinds_x = "stable"
      kind_w = "stable"
      XdX = X
      
      call pt_F_NP(model, z, np, 0.0_pr, kinds_x, kind_w, XdX, ns, S, F, df, Vl=Vl, Vw=Vw)
      ! print *, F
      ! print "(*(I10,2x))", (i, i=1, psize)
      do i=1,nc*np+np + 2
         XdX = X
         dX = XdX(i) * eps
         XdX(i) = XdX(i) + dX

         call pt_F_NP(model, z, np, 0.0_pr, kinds_x, kind_w, XdX, ns, S, F1, tmp, Vl=Vl, Vw=Vw)
         XdX(i) = XdX(i) - 2*dx
         call pt_F_NP(model, z, np, 0.0_pr, kinds_x, kind_w, XdX, ns, S, F2, tmp, Vl=Vl, Vw=Vw)
         dfnum(:, i) = (F1 - F2)/(2*dX)

         ! print *, i
         ! print fmt, df(:, i)
         ! print fmt, dfnum(:, i)
         ! print *, ""

         ! if (i > nc*np) print *, "=============================================="
      end do

      loc = maxloc(abs(df - dfnum))
      ! print *, maxval(abs((df - dfnum))), maxloc(abs(df - dfnum)), df(loc(1), loc(2)), dfnum(loc(1), loc(2))
      maxerr = maxval(abs((df - dfnum)))


   end subroutine numdiff
end program main
