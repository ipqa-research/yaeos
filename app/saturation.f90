program main
   use yaeos
   use fortime, only: timer
   integer, parameter :: nc=2
   type(NRTL) :: ge_model

   real(pr), dimension(nc,nc) :: a, b, c
   real(pr), dimension(nc) :: tc, pc, w
   type(MHV) :: mixrule
   type(CubicEOS) :: eos
   type(EquilibriumState) :: eq
   type(PXEnvel2) :: px
   type(timer) :: t
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

   ! print *, "Pxy a mano"


   print *, "Pxy envlop"
   n = [0.001_pr, 0.999_pr]
   eq = saturation_pressure(eos, n, 273._pr + 100._pr, kind="bubble")
   if (eq%iters > 1000) error stop 1


   n = [0, 1]
   call t%timer_start()
   px = px_envelope_2ph(&
      eos, z0=n, alpha0=0.001_pr, &
      z_injection=[1.0_pr, 0.0_pr], first_point=eq, iterations=30, delta_0=0.2_pr &
   )
   call t%timer_stop()

   do i=1,size(px%points)
      write(2, *) px%points(i)%x(1), px%points(i)%y(1), px%points(i)%P
   end do

   print *, px%cps
end program main