program main
   use yaeos, only: pr, CubicEoS
   use fixtures_models, only: binary_NRTL_SRK_HV
   use auxiliar_functions, only: allclose
   use testing_aux, only: test_ok, test_title

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
   call assert(allclose(dDij, test_dDij, test_tol), &
      "Huron-Vidal Mixing Rule: dDij does not match expected value")

   ! Below are the numerical derivatives used to generate the test values 

   
   ! print *, "numdiff"

   ! print *, "dDdT", dDdT,   (Dhv(n, T+0.01) - Dhv(n, T))/0.01
   ! print *, "dDdT2", dDdT2, (Dhv(n, T+0.01) - 2*Dhv(n, T) + Dhv(n, T-0.01))/(0.01**2)

   ! print *, "dDi"
   ! do i = 1, nc
   !   dn = 0
   !   dn(i) = 0.001
   !   print *, dDi(i), (Dhv(n + dn, T) - Dhv(n - dn, T))/(2 * dn(i))
   ! end do

   ! print *, "dDidT"
   ! do i = 1, nc
   !    dn = 0
   !    dn(i) = 0.005

   !    print *, dDidT(i), (&
   !         Dhv(n + dn, T + dn(i)) &
   !       - Dhv(n + dn, T - dn(i)) &
   !       - Dhv(n - dn, T + dn(i)) &
   !       + Dhv(n - dn, T - dn(i)) &
   !       )/(4 * dn(i)**2)
   ! end do

   ! print *, "dDij"

   ! dn = 0
   ! dn(1) = 1e-5
   ! print *, dDij(1, 1), (Dhv(n + dn, T) - 2*Dhv(n, T) + Dhv(n - dn, T))/(dn(1)**2)

   ! dn = 0
   ! dn(2) = 1e-5
   ! print *, dDij(2, 2), (Dhv(n + dn, T) - 2*Dhv(n, T) + Dhv(n - dn, T))/(dn(2)**2)
  
   ! dx = 1e-6
   
   ! print *, dDij(1, 2), (&
   !         Dhv(n + [dx, dx], T) &
   !       - Dhv(n + [-dx, dx], T) &
   !       - Dhv(n + [dx, -dx], T) &
   !       + Dhv(n + [-dx, -dx], T) &
   !       )/(4 * dx **2)
  

   ! contains

   !real(pr) function Dhv(n, T)
   !   real(pr), intent(in) :: n(:), T

   !   real(pr) :: Tc(size(n)), Tr(size(n))
   !   real(pr), dimension(size(n)) :: ai, daidt, dadit2
   !   real(pr) :: D, dDdT, dDdT2, dDi(size(n)), dDidT(size(n)), dDij(size(n), size(n))
   !   
   !   Tc = model%components%Tc
   !   Tr = T/Tc

   !   call model%alpha%alpha(Tr, ai, daidt, daidt2)

   !   ai = ai*model%ac
   !   daidt = daidt*model%ac/Tc
   !   daidt2 = daidt2*model%ac/Tc**2

   !   call model%mixrule%Dmix(n, T, ai, daidt, daidt2, D, dDdT, dDdT2, dDi, dDidT, dDij)
   !   Dhv = D
   !end function
end program main
