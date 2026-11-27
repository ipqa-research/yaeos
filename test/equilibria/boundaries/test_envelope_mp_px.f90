program main
   !! Test for multi-phase envelope calculation
   !! In this test we calculate two phase envelope from a previously known
   !! double saturation point.
   use testing_aux, only: assert, test_title
   use yaeos__constants, only: pr
   use yaeos
   implicit none
   integer, parameter :: nc = 15, np=3
   integer, parameter :: psize = np*nc + np + 2

   real(pr) :: z(nc)

   real(pr) :: P, T, X(psize), betas(np)
   real(pr) :: dF(psize, psize), F(psize), dXdS(psize)

   real(pr) :: x_l(np, nc), w(nc), tmp(nc)
   real(pr) :: K(np, nc)
   real(pr) :: beta_w
   real(pr) :: z0(nc), zi(nc)

   integer :: ns
   real(pr) :: S

   integer :: i, j, lb, ub

   integer :: iters

   type(CubicEoS) :: model
   type(PTEnvelMP) :: pt
   type(PXEnvelMP) :: px
   character(len=14) :: kinds_x(np)
   character(len=14) :: kind_w
   real(pr) :: Vl(np), Vw


   print *, test_title("Multi-phase PX envelope test")
   kinds_x = "stable"
   kind_w = "stable"

   model = get_model()

   ! Seed near the known four-phase point: the main phase x1 = z, a
   ! methane-rich liquid x2 (beta = 0, incipient), a heavy liquid x3 rich in
   ! the last pseudo-component, and the incipient vapor w. The values are the
   ! converged solution rounded to four significant digits. A cruder guess
   ! (x3 = 1e-10 for every light component, P = 70 bar) made Newton stall
   ! away from the solution on some toolchains (Ubuntu 26.04 amd64v3 builds).
   x_l(1, :) = z
   x_l(2, :) = [2.556e-02, 8.141e-03, 8.748e-01, 6.529e-02, 1.731e-02, &
      1.289e-03, 4.599e-03, 9.125e-04, 1.297e-03, 8.108e-04, &
      2.873e-05, 1.174e-11, 4.841e-13, 4.636e-06, 1.915e-23]
   x_l(3, :) = [1.081e-03, 6.104e-03, 1.492e-01, 4.668e-02, 2.442e-02, &
      2.121e-03, 1.216e-02, 3.096e-03, 5.719e-03, 4.608e-03, &
      2.610e-02, 2.374e-03, 1.101e-07, 1.374e-15, 7.163e-01]
   w = [1.007e-02, 9.510e-03, 6.284e-01, 1.082e-01, 5.291e-02, &
      6.007e-03, 2.550e-02, 7.919e-03, 1.297e-02, 1.421e-02, &
      3.512e-02, 2.299e-04, 7.104e-05, 8.890e-02, 1.979e-10]

   z0 = z
   zi = 0
   zi(2) = 1

   ! Initial guess for betas, P and T
   betas = [0.9989, 0.0, 1.054e-3]
   P = 108
   T = 261.8

   ns = psize

   X(1:nc) = log(x_l(1,:)/w)
   X(nc+1:2*nc) = log(x_l(2,:)/w)
   X(2*nc+1:3*nc) = log(x_l(3,:)/w)
   X(3*nc+1:3*nc+np) = betas
   X(3*nc + np + 1) = log(P)
   X(3*nc + np + 2) = log(T)

   ! Specify beta of the second main phase: it is incipient (beta = 0)
   ns = 3*nc+2
   S = 0

   ! Find an initial point by solving a known 4ph point
   call pt%solve_point(&
      model, z, np=np, beta_w=0.0_pr, kinds_x=kinds_x, kind_w=kind_w, &
      X=X, ns=ns, S=S, dXdS=dXdS, &
      F=F, dF=dF, Vl=Vl, Vw=Vw, iters=iters, max_iterations=100)
   P = exp(X(3*nc+np+1))
   T = exp(X(3*nc+np+2))
   betas = X(3*nc+1:3*nc+np)

   ! alpha
   X(3*nc + np + 2) = 0.01
   ns = size(X)

   call numdiff

   px = px_envelope(&
      model, z0, zi, np, T=T, kinds_x=kinds_x, kind_w=kind_w, &
      x_l0=x_l, w0=w, betas0=betas, p0=P, alpha0=0.0_pr,&
      ns0=psize, dS0=1e-2_pr, beta_w=0.0_pr, points=200)
   call assert(maxval(abs(px%points(1)%betas - [0.99, 0.0, 1.05e-3])) < 1e-2, "First point betas")
   call assert(abs(px%points(1)%P - 108.015 )< 1e-2, "First point P")
   call assert(abs(px%points(1)%T - 261.828 )< 1e-2, "First point T")
   ! Every point of the line is a true equilibrium state
   do j = 1, size(px%points)
      associate(pt => px%points(j))
         call check_equilibrium(pt%x_l, pt%w, pt%betas, pt%P, pt%T, pt%kinds_x, pt%kind_w, &
            px%alpha(j)*zi + (1 - px%alpha(j))*z0, "Envelope point")
      end associate
   end do

   i = size(px%points)
   call assert(abs(px%points(i)%P) > 9, "End at low pressure")

contains

   subroutine check_equilibrium(x_l, w, betas, P, T, kinds_x, kind_w, zref, name)
      !! Independent check of a converged point: isofugacity of every main
      !! phase with the incipient phase, and the material balance.
      real(pr), intent(in) :: x_l(np, nc) !! Main phases compositions
      real(pr), intent(in) :: w(nc) !! Incipient phase composition
      real(pr), intent(in) :: betas(np) !! Main phases fractions
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]
      character(len=14), intent(in) :: kinds_x(np) !! Main phases kinds
      character(len=14), intent(in) :: kind_w !! Incipient phase kind
      real(pr), intent(in) :: zref(nc) !! Global composition of the point
      character(len=*), intent(in) :: name
      real(pr) :: lnphi_x(nc)
      real(pr) :: lnphi_w(nc)
      real(pr) :: zcalc(nc)
      integer :: l

      call model%lnphi_pt(w, P, T, root_type=kind_w, lnPhi=lnphi_w)
      zcalc = 0
      do l=1,np
         call model%lnphi_pt(x_l(l, :), P, T, root_type=kinds_x(l), lnPhi=lnphi_x)
         call assert(maxval(abs(log(x_l(l, :)) + lnphi_x - log(w) - lnphi_w)) < 1e-6_pr, &
            name // ": isofugacity of every main phase with the incipient phase")
         zcalc = zcalc + betas(l)*x_l(l, :)
      end do
      call assert(maxval(abs(zcalc - zref)) < 1e-8_pr, name // ": material balance")
   end subroutine check_equilibrium

   subroutine numdiff
      use yaeos__equilibria_boundaries_phase_envelopes_mp_px, only: px_F_NP
      real(pr) :: dfnum(psize, psize), F1(psize), F2(psize)
      real(pr) :: tmp(psize, psize)
      real(pr) :: XdX(psize)
      real(pr) :: eps = 1e-3, dx
      character(len=*), parameter :: fmt="(*(E10.2,2x))"

      integer :: i, j, loc(2)

      call px_F_NP(model, z0, zi, np, T, 0.0_pr, kinds_x, kind_w, X, ns, S, F, dF)

      XdX = X
      do i=1,nc*np+np + 2
         XdX = X
         dX = maxval([abs(XdX(i) * eps), 1e-5_pr])
         XdX(i) = XdX(i) + dX
         call px_F_NP(model, z0, zi, np, T, 0.0_pr, kinds_x, kind_w, Xdx, ns, S, F1, tmp)
         XdX(i) = XdX(i) - 2*dx
         call px_F_NP(model, z0, zi, np, T, 0.0_pr, kinds_x, kind_w, Xdx, ns, S, F2, tmp)
         dfnum(:, i) = (F1 - F2)/(2*dX)

      end do

      where (dfnum == 0)
         dfnum = 1e-15_pr
      end where

      where (df == 0)
         df = 1e-15_pr
      end where

      call assert(maxval(abs((df - dfnum))) < 2e-1, "Numerical derivative matches analytical")
   end subroutine numdiff

   type(CubicEoS) function get_model()
      real(pr) :: tc(15), pc(15), w(15), kij(15, 15)
      z = [0.0048, 0.00919, 0.43391, 0.1101, 0.06544, 0.00789, 0.03787, &
         0.01279, 0.02248, 0.02698, 0.22738, 0.03747, 0.0023, 0.00054, 0.00086]
      Tc = [126.2, 304.2, 190.6, 305.4, 369.8, 408.1, 425.2, 460.4, &
         469.6, 507.4, 691.81, 956.77, 1118.6, 1325.03, 1445.73]
      Pc = [33.94, 73.76, 46., 48.84, 42.46, 36.48, 38., 33.84, 33.74, 29.69, 19.46, 13.08, 10.66, 10.28, 17.3]
      w = [0.04, 0.225, 0.008, 0.098, 0.152, 0.176, 0.193, 0.227, 0.251, 0.296, 0.68, 1.208, 0.949, 0.182, 1.274]

      kij(1, :) = [0., -0.032, 0.028, 0.041, 0.076, 0.094, 0.07, 0.087, 0.088, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08]
      kij(2, :) = [-0.032, 0., 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.1, 0.1, 0.1, 0.1, 0.1]
      kij(3, :) = [0.028, 0.12, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.017]
      kij(4, :) = [0.041, 0.12, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.017]
      kij(5, :) = [0.076, 0.12, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.017]
      kij(6, :) = [0.094, 0.12, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.017]
      kij(7, :) = [0.07, 0.12, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.017]
      kij(8, :) = [0.087, 0.12, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.017]
      kij(9, :) = [0.088, 0.12, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.017]
      kij(10, :) = [0.08, 0.12, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.017]
      kij(11, :) = [0.08, 0.1, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
      kij(12, :) = [0.08, 0.1, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
      kij(13, :) = [0.08, 0.1, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
      kij(14, :) = [0.08, 0.1, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
      kij(15, :) = [0.08, 0.1, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0., 0., 0., 0., 0.]

      get_model = PengRobinson76(Tc, Pc, w, kij=kij)
   end function get_model
end program main
