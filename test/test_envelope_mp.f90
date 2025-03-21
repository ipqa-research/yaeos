program main
   !! Test for multi-phase envelope calculation
   !! In this test we calculate two phase envelope from a previously known
   !! double saturation point.
   use testing_aux, only: assert, test_title
   use yaeos
   use yaeos__equilibria_boundaries_phase_envelopes_mp, only: pt_F_NP, solve_point, pt_envelope, PTEnvelMP
   implicit none
   integer, parameter :: nc = 15, np=3
   integer, parameter :: psize = np*nc + np + 2

   real(pr) :: z(nc)

   real(pr) :: P, T, X(psize), betas(np)
   real(pr) :: dF(psize, psize), F(psize)

   real(pr) :: x_l(np, nc), w(nc), tmp(nc)
   real(pr) :: K(np, nc)

   integer :: ns
   real(pr) :: S

   integer :: i, lb, ub

   integer :: iters

   type(CubicEoS) :: model
   type(PTEnvelMP) :: env


   print *, test_title("Multi-phase envelope test")

   model = get_model()

   ! Build composition for three main phases and incipient one
   x_l(1, :) = z
   x_l(2, :) = [2.26748304e-02, 8.77274788e-03, 8.56790518e-01, 7.53110308e-02, 2.25791412e-02, 1.83124730e-03, 6.75557205e-03, 1.46973097e-03, 2.15197976e-03, 1.51272521e-03, 1.15635706e-04, 8.60487017e-11, 8.01982281e-12, 3.48403955e-05, 3.09839163e-21]
   x_l(3, :) = 1e-10
   x_l(3, nc) = 1 - 1e-10*(nc-1)
   w = [1.00840444e-02, 9.76751155e-03, 6.35587304e-01, 1.11530262e-01, 5.50813917e-02, 6.27551663e-03, 2.67981694e-02, 8.36513089e-03, 1.37645767e-02, 1.51732534e-02, 3.84679648e-02, 1.63919981e-04, 5.56318734e-05, 6.88853217e-02, 8.30156039e-10]

   ! Initial guess for betas
   betas = [0.98, 0.01, 0.]
   betas = betas/sum(betas)


   ! Build variables vector
   do i=1, np
      lb = (i-1)*nc + 1
      ub = i*nc
      X(lb:ub) = log(x_l(i, :)/w)
   end do
   X(3*nc+1:np*nc+np) = betas
   P = 70
   T = 260
   X(nc*np + np + 2) = log(T)
   X(nc*np + np + 1) = log(P)

   ! Select specification variable
   ns = nc*np + np + 1
   S = X(ns)

   betas = X(3*nc+1:np*nc+np)

   betas = [betas(1), 0._pr, betas(3)]
   env = pt_envelope(model, z, np, x_l, w, betas, P, T, np*nc+2, 0.01_pr)

   call assert(maxval(abs(env%points(1)%betas - [0.99, 0.0, 1.05e-3])) < 1e-2, "First point betas")
   call assert(abs(env%points(1)%P - 108.015 )< 1e-2, "First point P")
   call assert(abs(env%points(1)%T - 261.828 )< 1e-2, "First point T")

   i = size(env%points)
   call assert(abs(env%points(i)%P) < 1, "End at low pressure")

contains

   type(CubicEoS) function get_model()
      real(pr) :: tc(15), pc(15), w(15), kij(15, 15)
      z = [0.0048, 0.00919, 0.43391, 0.1101, 0.06544, 0.00789, 0.03787, 0.01279, 0.02248, 0.02698, 0.22738, 0.03747, 0.0023, 0.00054, 0.00086]
      Tc = [126.2, 304.2, 190.6, 305.4, 369.8, 408.1, 425.2, 460.4, 469.6, 507.4, 691.81, 956.77, 1118.6, 1325.03, 1445.73]
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
