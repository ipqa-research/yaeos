program main
   use yaeos
   use testing_aux, only: assert, test_title
   use yaeos__math, only: intersect_one_line, point

   integer, parameter :: nc=16
   integer, parameter :: np=1

   type(CubicEoS) :: ar_model

   real(pr) :: tc(nc), pc(nc), w(nc)

   type(QMR) :: mixrule
   real(pr) :: kij(nc, nc), lij(nc, nc)

   real(pr) :: z(nc)

   type(PTEnvelMP) :: env

   real(pr) :: x_l0(np, nc), w0(nc), betas0(np), p0, t0
   character(len=14) :: kinds_x(np), kind_w

   type(Point), allocatable :: intersection(:)


   type(EquilibriumState) :: sat

   write(*, *) test_title("PT with self crossing")


   tc = [304.2_pr, 126.2_pr, 190.6_pr, 305.4_pr, 369.8_pr, 408.1_pr, 425.2_pr, 460.4_pr, 469.6_pr, 506.35_pr, 566.55_pr, 647.06_pr, 719.44_pr, 784.93_pr, 846.33_pr, 919.39_pr]
   pc = [72.8_pr, 33.5_pr, 45.4_pr, 48.2_pr, 41.9_pr, 36.0_pr, 37.5_pr, 33.4_pr, 33.3_pr, 33.9_pr, 25.3_pr, 19.1_pr, 14.2_pr, 10.5_pr, 7.5_pr, 4.76_pr]
   w = [0.225_pr, 0.04_pr, 0.008_pr, 0.098_pr, 0.152_pr, 0.176_pr, 0.193_pr, 0.227_pr, 0.251_pr, 0.299_pr, 0.3884_pr, 0.5289_pr, 0.6911_pr, 0.8782_pr, 1.1009_pr, 1.4478_pr]
   z = [0.710319, 0.001392, 0.04727, 0.011687, 0.008613, 0.001044, 0.009541, 0.004582, 0.006235, 0.009628, 0.05340408, 0.04752839, 0.03690424, 0.02809752, 0.01705432, 0.00670045]

   ar_model = PengRobinson78(tc, pc, w)

   kij(1, :) = [0.0_pr, -0.02_pr, 0.075_pr, 0.08_pr, 0.08_pr, 0.085_pr, 0.085_pr, 0.085_pr, 0.085_pr, 0.095_pr, 0.095_pr, 0.095_pr, 0.095_pr, 0.095_pr, 0.095_pr, 0.095_pr]
   kij(2, :) = [-0.02_pr, 0.0_pr, 0.08_pr, 0.07_pr, 0.07_pr, 0.06_pr, 0.06_pr, 0.06_pr, 0.06_pr, 0.05_pr, 0.1_pr, 0.12_pr, 0.12_pr, 0.12_pr, 0.12_pr, 0.12_pr]
   kij(3, :) = [0.075_pr, 0.08_pr, 0.0_pr, 0.003_pr, 0.01_pr, 0.018_pr, 0.018_pr, 0.025_pr, 0.026_pr, 0.036_pr, 0.049_pr, 0.073_pr, 0.098_pr, 0.124_pr, 0.149_pr, 0.181_pr]
   kij(4, :) = [0.08_pr, 0.07_pr, 0.003_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   kij(5, :) = [0.08_pr, 0.07_pr, 0.01_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   kij(6, :) = [0.085_pr, 0.06_pr, 0.018_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   kij(7, :) = [0.085_pr, 0.06_pr, 0.018_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   kij(8, :) = [0.085_pr, 0.06_pr, 0.025_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   kij(9, :) = [0.085_pr, 0.06_pr, 0.026_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   kij(10, :) = [0.095_pr, 0.05_pr, 0.036_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   kij(11, :) = [0.095_pr, 0.1_pr, 0.049_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   kij(12, :) = [0.095_pr, 0.12_pr, 0.073_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   kij(13, :) = [0.095_pr, 0.12_pr, 0.098_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   kij(14, :) = [0.095_pr, 0.12_pr, 0.124_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   kij(15, :) = [0.095_pr, 0.12_pr, 0.149_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   kij(16, :) = [0.095_pr, 0.12_pr, 0.181_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]

   lij(1, :) = [0.0_pr, -0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   lij(2, :) = [-0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   lij(3, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   lij(4, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   lij(5, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   lij(6, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   lij(7, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   lij(8, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   lij(9, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   lij(10, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   lij(11, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   lij(12, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   lij(13, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   lij(14, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   lij(15, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
   lij(16, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]

   mixrule = QMR(k=kij, l=lij)


   call ar_model%set_mixrule(mixrule)


   sat = saturation_temperature(ar_model, z, P=0.05_pr, kind="dew", t0=300._pr)

   x_l0(1, :) = z
   w0 = sat%x
   betas0 = 1
   kinds_x = "vapor"
   kind_w = "liquid"
   p0 = sat%P
   t0 = sat%T

   env = pt_envelope(& 
      ar_model, z, np, kinds_x, kind_w, x_l0, w0, betas0, p0, t0, &
      ns0=nc*np+np+1, ds0=1e-1_pr, beta_w=0._pr &
   )

   intersection = intersect_one_line(env%points%T, env%points%P)

   call assert(size(env%Pc) == 2, "Two critical points found")
   call assert(env%points(size(env%points))%P > 1500._pr, "Envelope should end at high pressure")
   call assert(size(intersection) == 1, "Envelope should have one self intersection")
   call assert(abs(intersection(1)%x - 305._pr) < 1._pr, "Intersection should be near 305 K")
   call assert(abs(intersection(1)%y - 330._pr) < 1._pr, "Intersection should be near 330 bar")

end program