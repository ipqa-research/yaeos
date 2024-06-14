program main
   !! Example of using CubicEoS with Huron-Vidal mixing rules with an
   !! NRTL model as the Ge model

   use yaeos, only: pr, SoaveRedlichKwong, CubicEoS, NRTL
   use yaeos__models_cubic_mixing_rules_huron_vidal, only: HV
   
   implicit none
   integer, parameter :: nc = 2

   real(pr) :: T=150, n(nc)
   
   real(pr) :: a(nc, nc), b(nc, nc), c(nc, nc) ! NRTL parameters
   real(pr) :: tc(nc), pc(nc), w(nc) ! Cubic EoS parameters

   real(pr) :: alpha(nc), Tr(nc)
   type(NRTL) :: ge_model ! Excess Gibbs model that will be used
   type(CubicEoS) :: model ! Main model
   type(HV) :: mixrule

   real(pr), dimension(nc) :: ai, daidt, daidt2
   real(pr) ::  D, dDdT, dDdT2, dDi(nc), dDidT(nc), dDij(nc,nc)

   tc = [647.13, 514.0]
   w = [0.344861, 0.643558]
   pc = [220.55, 61.37]
   a = 0; b = 0; c = 0

   a(1, 2) = 3.458
   a(2, 1) = -0.801

   b(1, 2) = -586.1
   b(2, 1) = 246.2

   c(1, 2) = 0.3
   c(2, 1) = 0.3
   ge_model = NRTL(a, b, c)

   model = SoaveRedlichKwong(tc, pc, w)
   mixrule = HV(ge_model, model%b)
   mixrule%q = 0.53

   Tr = T/Tc
   call model%alpha%alpha(Tr, alpha, daidt, daidt2)

   n = [0.2, 0.8]

   ai = alpha * model%ac
   daidt = alpha * model%ac / Tc
   daidt2 = alpha * model%ac / Tc**2
   call mixrule%Dmix(n, T, ai, daidt, daidt2, D, dDdT, dDdT2, dDi, dDidT, dDij)
   print *, D
end program
