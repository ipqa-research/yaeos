program main
   use yaeos, only: pr, CubicEoS
   use fixtures_models, only: binary_NRTL_SRK_HV
   use auxiliar_functions, only: allclose
   use testing_aux, only: test_title, assert

   integer, parameter :: nc = 2
   real(pr) :: test_tol = 1e-7
   real(pr) :: test_D =    30.622979288145512
   real(pr) :: test_dDdT = -9.2515493337243210E-002
   real(pr) :: test_dDdT2 =4.8515216513941623E-004
   real(pr) :: test_dDi(nc) = [44.282156269357145, 65.486908012229634]
   real(pr) :: test_dDidT(nc) = [-0.11932459578485301, -0.20145758095042413]
   real(pr) :: test_dDij(nc, nc) = reshape( &
      [26.205198410827229, 48.801394909170199, &
      48.801394909170199,69.658285068205771 ], [nc, nc])

   real(pr) :: ai(nc), daidt(nc), daidt2(nc)
   real(pr) :: n(nc), T, Tr(nc), Tc(nc), dn(nc)

   real(pr) :: D, dDdT, dDdT2, dDi(nc), dDidT(nc), dDij(nc,nc), dx

   type(CubicEoS) :: model
   integer :: i, j

   model = binary_NRTL_SRK_HV()

   write(*, *) test_title("HURON_VIDAL MIXING RULE")

   n = [0.2, 0.8]
   T = 150
   Tc = model%components%Tc
   Tr = T/Tc

   call model%alpha%alpha(Tr, ai, daidt, daidt2)

   ai = ai*model%ac
   daidt = daidt*model%ac/Tc
   daidt2 = daidt2*model%ac/Tc**2

   call model%mixrule%Dmix(n, T, ai, daidt, daidt2, D, dDdT, dDdT2, dDi, dDidT, dDij)

   call assert(allclose([D], [test_D], test_tol), &
      "Huron-Vidal Mixing Rule: D does not match expected value")

   call assert(allclose([dDdT], [test_dDdT], test_tol), &
      "Huron-Vidal Mixing Rule: dDdT does not match expected value")
   call assert(allclose([dDdT2], [test_dDdT2], test_tol), &
      "Huron-Vidal Mixing Rule: dDdT2 does not match expected value")
   call assert(allclose(dDi, test_dDi, test_tol), &
      "Huron-Vidal Mixing Rule: dDi does not match expected value")
   call assert(allclose(dDidT, test_dDidT, test_tol), &
      "Huron-Vidal Mixing Rule: dDidT does not match expected value")
   call assert(allclose([dDij], [test_dDij], test_tol), &
      "Huron-Vidal Mixing Rule: dDij does not match expected value")

   call test_kij

contains
   subroutine test_kij
      use yaeos, only: pr, CubicEoS, NRTLHV, init_hvnrtl, HV_NRTL, PengRobinson76
      integer, parameter :: nc = 3
      type(CubicEoS) :: model_kij
      type(HV_NRTL) :: mixrule
      type(NRTLHV) :: ge
      real(pr) :: kij(nc, nc)
      logical :: use_kij(nc, nc)
      real(pr) :: alpha(nc, nc), gji(nc, nc)
      real(pr) :: Tc(nc), Pc(nc), w(nc)

      real(pr) :: D, dDi(nc), dDT, dDdT2, dDidT(nc), dDij(nc, nc)
      real(pr) :: Ar1, Ar2
      real(pr) :: n(nc), V, T

      n = [2, 3, 5]
      V = 1
      T = 300

      use_kij = .false.
      use_kij(1, 2) = .true.
      use_kij(2, 1) = .true.

      alpha(1, :) = [0.0, 0.2, 0.3]
      alpha(2, :) = [0.3, 0.0, 0.1]
      alpha(3, :) = [0.5, 0.0, 0.1]

      gji(1, :) = [0.0, 0.1, 0.2]*100
      gji(2, :) = [0.2, 0.0, 0.3]*100
      gji(3, :) = [0.2, 0.9, 0.0]*100

      kij = 0.0
      kij(1, 2) = 0.1
      kij(2, 1) = 0.1

      ! CO2, Methane, Butane
      Tc =  [304.21, 190.564, 425.12]
      Pc =  [73.83000000000001, 45.99, 37.96]
      w =  [0.223621, 0.0115478, 0.200164]

      model_kij = PengRobinson76(Tc, Pc, w, kij=kij)
      use_kij = .true.

      use_kij(1, 2) = .true.
      use_kij(2, 1) = .true.

      ge = NRTLHV(b=model_kij%b, alpha=alpha, gji0=gji, gjiT=0*gji)
      mixrule = init_hvnrtl(b=model_kij%b, del1=model_kij%del1, alpha=alpha, gji0=gji, gjiT=0*gji, use_kij=use_kij, kij=kij)

      call model_kij%residual_helmholtz(n, V, T, Ar=Ar1)

      call model_kij%set_mixrule(mixrule)
      call model_kij%residual_helmholtz(n, V, T, Ar=Ar2)

      call assert(abs(Ar1- Ar2) < 1e-10, &
         "Huron-Vidal mixing rule: Ar must match with the one calculated with QMR")
   end subroutine test_kij
end program main
