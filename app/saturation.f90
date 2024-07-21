program main
   use yaeos
   integer, parameter :: nc=2
   type(NRTL) :: ge_model

   real(pr), dimension(nc,nc) :: a, b, c
   real(pr), dimension(nc) :: tc, pc, w
   type(MHV) :: mixrule
   type(CubicEOS) :: eos
   type(EquilibriumState) :: eq
   real(pr) :: n(nc)

   integer :: i

   tc = [647.13999999999999, 513.91999999999996]
   pc = [220.63999999999999, 61.479999999999997]
   w =  [0.34399999999999997, 0.64900000000000002]
   a = 0; b = 0; c = 0

   ! NRTL model parameters
   a(1, 2) = 3.458
   a(2, 1) = -0.801

   b(1, 2) = -586.1
   b(2, 1) = 246.2

   c(1, 2) = 0.3
   c(2, 1) = 0.3

   eos = PengRobinson76(tc, pc, w)
   ge_model = NRTL(a, b, c)
   mixrule = MHV(ge=ge_model, q=-0.53_pr, b=eos%b)

   deallocate(eos%mixrule)
   eos%mixrule = mixrule

   do i=1,99
    n(1) = real(i, pr)/100
    n(2) = 1 - n(1)
    eq = saturation_pressure(eos, n, 273._pr + 100._pr, kind="bubble")
    print *, eq%x(1), eq%y(1), eq%P
   end do

   call ge_model%ln_activity_coefficient([0.3_pr, 0.7_pr], 250._pr, n)
   print *, n
end program main


! for i, T in enumerate(np.linspace(50+273, 200+273, 5)):
!     i=4
!     xs = np.linspace(0.001, 0.999, 100)
!     ys = []
!     ps = []
!
!     for x1 in xs:
!         x = [x1, 1-x1]
!         p, x, y, vx, vy, beta = yaeos.yaeos_c.saturation_pressure(model.id, x, T, "bubble")
!         ps.append(p)
!         ys.append(y[0])
!
!     plt.plot(xs, ps, color=colors[i])
!     plt.plot(ys, ps, color=colors[i])
! end program
!
