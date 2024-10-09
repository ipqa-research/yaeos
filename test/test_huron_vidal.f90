program main
   use yaeos, only: pr, CubicEoS
   use fixtures_models, only: binary_NRTL_SRK_HV
   use auxiliar_functions, only: allclose

   integer, parameter :: nc = 2
   real(pr) :: test_tol = 1e-7
   real(pr) :: test_D =    30.622979288145512
   real(pr) :: test_dDdT = -9.2515493337243210E-002
   real(pr) :: test_dDdT2 =4.8515216513941623E-004
   real(pr) :: test_dDi(nc) = [44.282156269357145, 65.486908012229634]
   real(pr) :: test_dDidT(nc) = [-0.11932459578485301, -0.20145758095042413]
   real(pr) :: test_dDij(nc, nc) = reshape( &
      [25.301358190454181, 49.027354964263466, &
      49.027354964263466, 69.601795054432472], [nc, nc])

   real(pr) :: ai(nc), daidt(nc), daidt2(nc)
   real(pr) :: n(nc), T, Tr(nc), Tc(nc)

   real(pr) :: D, dDdT, dDdT2, dDi(nc), dDidT(nc), dDij(nc,nc)

   type(CubicEoS) :: model

   model = binary_NRTL_SRK_HV()

   print *, "HURON_VIDAL MIXING RULE"

   n = [0.2, 0.8]
   T = 150
   Tc = model%components%Tc
   Tr = T/Tc

   call model%alpha%alpha(Tr, ai, daidt, daidt2)

   ai = ai*model%ac
   daidt = daidt*model%ac/Tc
   daidt2 = daidt2*model%ac/Tc**2

   call model%mixrule%Dmix(n, T, ai, daidt, daidt2, D, dDdT, dDdT2, dDi, dDidT, dDij)
   print *, D
   print *, dDdT
   print *, dDdT2
   print *, dDi
   print *, dDidT
   print *, dDij

   if (.not. allclose([D], [test_D], test_tol)) error stop 1
   if (.not. allclose([dDdT], [test_dDdT], test_tol)) error stop 1
   if (.not. allclose([dDdT2], [test_dDdT2], test_tol)) error stop 1
   if (.not. allclose([dDi], [test_dDi], test_tol)) error stop 1
   if (.not. allclose([dDidT], [test_dDidT], test_tol)) error stop 1
   if (.not. allclose([dDij], [test_dDij], test_tol)) error stop 1



end program main
