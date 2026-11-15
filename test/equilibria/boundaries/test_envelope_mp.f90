program main
   !! Test for multi-phase envelope calculation
   !! In this test we calculate two phase envelope from a previously known
   !! double saturation point.
   use testing_aux, only: assert, test_title
   use yaeos
   use yaeos__equilibria_boundaries_phase_envelopes_mp, only: &
      pt_F_NP, pt_envelope, PTEnvelMP
   implicit none
   integer, parameter :: nc = 15, np=3
   integer, parameter :: psize = np*nc + np + 2

   real(pr) :: z(nc)

   real(pr) :: P, T, X(psize), betas(np)
   real(pr) :: dF(psize, psize), F(psize)

   real(pr) :: x_l(np, nc), w(nc), tmp(nc)
   real(pr) :: K(np, nc)
   real(pr) :: beta_w

   integer :: ns
   real(pr) :: S

   integer :: i, j, lb, ub

   integer :: iters

   type(CubicEoS) :: model
   type(PTEnvelMP) :: env
   character(len=14) :: kinds_x(np)
   character(len=14) :: kind_w


   print *, test_title("Multi-phase envelope test")

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

   ! Initial guess for betas, P and T
   betas = [0.9989, 0.0, 1.054e-3]
   P = 108
   T = 261.8
   kinds_x = ["liquid", "liquid", "liquid"]
   kind_w = "vapor"
   env = pt_envelope(model, z, np, kinds_x, kind_w, x_l, w, betas, P, T, np*nc+2, 0.001_pr, beta_w=0.0_pr)

   ! call env%write(1)

   call assert(maxval(abs(env%points(1)%betas - [0.99, 0.0, 1.05e-3])) < 1e-2, "First point betas")
   call assert(abs(env%points(1)%P - 108.015 )< 1e-2, "First point P")
   call assert(abs(env%points(1)%T - 261.828 )< 1e-2, "First point T")
   ! Every point of the line is a true equilibrium state
   do j = 1, size(env%points)
      associate(pt => env%points(j))
         call check_equilibrium(pt%x_l, pt%w, pt%betas, pt%P, pt%T, pt%kinds_x, pt%kind_w, &
            z, "Envelope point")
      end associate
   end do

   i = size(env%points)
   call assert(abs(env%points(i)%P) < 9, "End at low pressure")

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

   subroutine single
      ! Build variables vector
      do i=1, np
         lb = (i-1)*nc + 1
         ub = i*nc
         X(lb:ub) = log(x_l(i, :)/w)
      end do
      X(3*nc+1:np*nc+np) = betas
      X(nc*np + np + 2) = log(T)
      X(nc*np + np + 1) = log(P)

      ! Select specification variable
      ns = nc*np + np + 1
      S = X(ns)

      betas = X(3*nc+1:np*nc+np)

   end subroutine single
end program main
