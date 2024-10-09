program main
   use yaeos, only: pr, CubicEoS
   use fixtures_models, only: binary_NRTL_SRK_HV

   integer, parameter :: nc = 2
   real(pr) :: test_D = 30.797085158580309
   real(pr) :: test_dDdT = -8.4060598003984383E-002
   real(pr) :: test_dDdT2 = 3.5988912573413091E-004
   real(pr) :: test_dDi(nc) = [45.043635428115536, 65.731805235540278]
   real(pr) :: test_dDidT(nc) = [-8.2934765998547377E-002, -0.18941780037882283]
   real(pr) :: test_dDij(nc, nc) = reshape( &
      [25.390034959603035, 49.957034706240584, &
      49.957034706240584, 69.675496643514933], [nc, nc])

   real(pr) :: ai(nc), daidt(nc), daidt2(nc)
   real(pr) :: n(nc), T, Tr(nc), Tc(nc)

   real(pr) :: D, dDdT, dDdT2, dDi(nc), dDidT(nc), dDij(nc,nc)

   type(CubicEoS) :: model

   model = binary_NRTL_SRK_HV()

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
end program main
