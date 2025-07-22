program main
   use yaeos
   use yaeos__equilibria_multiphase_flash, only: pt_F_NP, solve_mp_flash_point&
      , pt_mp_flash, MPEquilibriumState
   use yaeos__equilibria_stability, only: min_tpd
   use yaeos__math, only: solve_system
   use testing_aux, only: test_title, assert
   use auxiliar_functions, only: allclose
   integer, parameter :: nc = 3, np = 2
   integer, parameter :: psize = nc * np + np + 3
   real(kind=pr) :: X(psize)
   real(kind=pr) :: z(nc), et, st
   type(CubicEoS) :: model
   type(EquilibriumState) :: fr

   integer :: ns1, ns2, i, iters, max_iters = 1000, beta_0_index
   logical :: less_phases
   real(kind=pr) :: S1, S2, tpd
   real(kind=pr) :: x_l(np, nc), w(nc), K(np, nc)
   real(kind=pr) :: betas(np), beta_w, P, T

   real(kind=pr) :: F(psize), dX(psize)
   real(kind=pr) :: df(psize, psize)

   character(len=14) :: kinds_x(np)
   character(len=14) :: kind_w

   integer :: j
   real(kind=pr) :: T0 = 150, Tf = 700
   real(kind=pr) :: P0 = 1, Pf = 350
   ! real(pr) :: T0=100, Tf=400
   ! real(pr) :: P0=5, Pf=400
   real(kind=pr) :: dT, dP
   integer, parameter :: npp = 500, nt = 500
   type(MPEquilibriumState) :: mpfr
   real(kind=pr) :: ts(nt, npp), ps(nt, npp), nps(nt, npp)

   dT = (Tf - T0) / nt
   dP = (Pf - P0) / npp

   write(*, *) test_title("MULTIPHASE FLASH SOLVER")

   model = get_model()

   T = 260
   P = 150
   w = [8.21852603e-01, 1.78147397e-01, 3.03035466e-16]
   x_l(1, :) = [2.64466950e-01, 1.33648382e-01, 6.01884668e-01]
   x_l(2, :) = [7.22697463e-02, 9.27730254e-01, 5.51527260e-13]
   beta_w = 0.01_pr
   betas = [0.66457915, 0.33542085]

   kinds_x = "liquid"
   kind_w = "vapor"

   ns1 = np * nc + np + 1 + 1
   S1 = log(P)
   ns2 = np * nc + np + 1 + 2
   S2 = log(T)
   
   block
      real(kind=pr) :: k(nc)
      k = 0.001_pr
      k(1) = 100
      fr = flash(model, z, T, P_spec=P, iters=i, k0=k)
   end block
   mpfr = pt_mp_flash(model, z, P, T)

   call assert(mpfr%np == 1, "Number of phases in multiphase flash")
   call assert(allclose(mpfr%x_l(1, :), fr%x, rtol=0.01_pr), "Heavy phase" // &
      " composition compared with SS flash")
   call assert(allclose(mpfr%w, fr%y, rtol=0.01_pr), "Light phase" // &
      " composition compared with SS flash")

   T = 260
   P = 50
   mpfr = pt_mp_flash(model, z, P, T)
   call assert(mpfr%np == 2, "Number of phases in multiphase flash")

   x_l(1, :) = [0.16813294321741362, 0.30875990831447953, 0.52310714846810680]
   x_l(2, :) = [0.11282300674090837, 0.88717246176596665, 4.5314931247530688E-006]
   w = [0.41554260326577303, 0.58445738903117439, 7.7030530774883127E-009]

   call assert(allclose(mpfr%betas, &
      [0.76466094081795022_pr, 8.7071006133493931E-002_pr, 0.14826805304855589_pr], &
      rtol=0.01_pr), "Beta")
   call assert(allclose(mpfr%x_l(1, :), x_l(1, :), rtol=0.01_pr), "x_l1" // &
      " phase composition")
   call assert(allclose(mpfr%x_l(2, :), x_l(2, :), rtol=0.01_pr), "x_l2" // &
      " phase composition")
   call assert(allclose(mpfr%w, w, rtol=0.01_pr), "w phase composition")
contains

   type(CubicEoS) function get_model()
      real(kind=pr) :: tc(nc), pc(nc), w(nc), kij(nc, nc)
      z = [0.2, 0.4, 0.4]
      tc = [190.564, 304.1282, 768.0]
      pc = [45.992, 73.773, 10.7]
      w = [0.01142, 0.22394, 0.8805]

      kij = 0.0_pr
      kij(1, 2) = 0.1
      kij(2, 1) = 0.1
      kij(2, 3) = 0.2
      kij(3, 2) = 0.2

      z = z/sum(z)

      get_model = PengRobinson78(Tc, Pc, w, kij=kij)
   end function get_model

   subroutine numdiff
      real(kind=pr) :: dfnum(psize, psize), F1(psize), F2(psize)
      real(kind=pr) :: tmp(psize, psize)
      real(kind=pr) :: XdX(psize)
      real(kind=pr) :: eps = 1e-2, dx
      character(len=*), parameter :: fmt = "(*(E11.3,2x))"
      integer :: i, j, loc(2)

      XdX = X
      print "(*(I10,2x))", (i, i = 1, psize)
      kinds_x = "stable"
      kind_w = "stable"
      do i = 1, nc * np + np + 3
         XdX = X
         dX = XdX(i) * eps
         XdX(i) = XdX(i) + dX

         call pt_F_NP(model, z, np, kinds_x, kind_w, Xdx, ns1, S1, ns2, S2, F1&
            , tmp)
         XdX(i) = XdX(i) - 2 * dx
         call pt_F_NP(model, z, np, kinds_x, kind_w, Xdx, ns1, S1, ns2, S2, F2&
            , tmp)
         dfnum(:, i) = (F1 - F2) / (2 * dX)

         ! print *, i
         ! print fmt, df(:, i)
         ! print fmt, dfnum(:, i)
         ! write(1, *) dfnum(:, i)
         ! print *, ""

         ! if (i > nc*np) print *,
         ! "=============================================="
      end do

      loc = maxloc(abs(df - dfnum))
      ! print *, maxval(abs((df - dfnum))), maxloc(abs(df - dfnum)), df(loc(1),
      ! loc(2)), dfnum(loc(1), loc(2))

   end subroutine numdiff

end program main
