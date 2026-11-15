program test_ge_excess_cp
   !! Excess heat capacity of `GeModel` and the Huron-Vidal NRTL
   !! temperature dependence `nrtl_hv_tdep`.
   !!
   !! - \(C_p^E = \partial H^E/\partial T\) checked against central finite
   !!   differences of `excess_enthalpy` (UNIFAC, a temperature dependent model).
   !! - \(C_p^E\) is consistent with \(-T\,\partial^2 G^E/\partial T^2\) and with
   !!   \(T\,\partial S^E/\partial T\).
   !! - \(\tau = g/(RT)\) and its derivatives; equivalence with the linear
   !!   dependence \(\tau = A + B/T\) for \(A=0,\ B=g/R\).
   use yaeos, only: pr, R, setup_unifac, UNIFAC, Groups
   use yaeos__models_ge_base, only: nrtl_hv_tdep, nrtl_hv_tdep_linear
   use testing_aux, only: assert, test_title
   implicit none

   integer, parameter :: nc = 3
   type(UNIFAC) :: model
   type(Groups) :: molecules(nc)
   real(pr) :: n(nc)
   real(pr) :: T
   real(pr) :: dT
   real(pr) :: Cpe
   real(pr) :: He_p
   real(pr) :: He_m
   real(pr) :: GeT2
   real(pr) :: SeT
   real(pr) :: Cpe_num

   real(pr) :: g(2, 2)
   real(pr) :: tau(2, 2)
   real(pr) :: dtau(2, 2)
   real(pr) :: dtau2(2, 2)
   real(pr) :: tau_l(2, 2)
   real(pr) :: dtau_l(2, 2)
   real(pr) :: dtau2_l(2, 2)
   real(pr) :: tau_p(2, 2)
   real(pr) :: tau_m(2, 2)
   real(pr) :: dtau_p(2, 2)
   real(pr) :: dtau_m(2, 2)

   write(*, *) test_title("GeModel excess Cp and NRTL-HV tau(T)")

   ! Ethane [CH3], Ethanol [CH3, CH2, OH], Methylamine [CH3NH2]
   molecules(1)%groups_ids = [1]
   molecules(1)%number_of_groups = [2]
   molecules(2)%groups_ids = [1, 2, 14]
   molecules(2)%number_of_groups = [1, 1, 1]
   molecules(3)%groups_ids = [28]
   molecules(3)%number_of_groups = [1]
   model = setup_unifac(molecules)

   n = [2.1_pr, 1.4_pr, 1.6_pr]
   T = 310._pr
   dT = 1e-3_pr

   call model%excess_Cp(n, T, Cpe)
   call model%excess_enthalpy(n, T + dT, He=He_p)
   call model%excess_enthalpy(n, T - dT, He=He_m)
   Cpe_num = (He_p - He_m)/(2*dT)
   call assert(abs(Cpe) > 0, "excess_Cp: UNIFAC mixture has a nonzero excess Cp")
   call assert(abs(Cpe - Cpe_num) < 1e-5_pr*abs(Cpe), "excess_Cp: Cpe equals dHe/dT (finite differences)")

   call model%excess_gibbs(n, T, GeT2=GeT2)
   call assert(abs(Cpe + T*GeT2) < 1e-10_pr*abs(Cpe), "excess_Cp: Cpe = -T d2Ge/dT2")

   call model%excess_entropy(n, T, SeT=SeT)
   call assert(abs(Cpe - T*SeT) < 1e-10_pr*abs(Cpe), "excess_Cp: Cpe = T dSe/dT")

   ! Ideal solution limit: a pure component has no excess Cp
   call model%excess_Cp([0._pr, 3._pr, 0._pr], T, Cpe)
   call assert(abs(Cpe) < 1e-10_pr, "excess_Cp: pure component has zero excess Cp")

   ! ==========================================================================
   ! nrtl_hv_tdep
   ! --------------------------------------------------------------------------
   g(1, :) = [0._pr, 1200._pr]
   g(2, :) = [-350._pr, 0._pr]

   call nrtl_hv_tdep(T, g, tau, dtau, dtau2)
   call assert(maxval(abs(tau - g/(R*T))) < 1e-13_pr*maxval(abs(tau)), "nrtl_hv_tdep: tau = g/(RT)")
   call assert(tau(1, 1) == 0 .and. tau(2, 2) == 0, "nrtl_hv_tdep: zero g gives zero tau")

   call nrtl_hv_tdep(T + dT, g, tau_p, dtau_p, dtau2)
   call nrtl_hv_tdep(T - dT, g, tau_m, dtau_m, dtau2)
   call nrtl_hv_tdep(T, g, tau, dtau, dtau2)
   call assert(maxval(abs((tau_p - tau_m)/(2*dT) - dtau)) < 1e-6_pr*maxval(abs(dtau)), &
      "nrtl_hv_tdep: dtau/dT vs finite differences")
   call assert(maxval(abs((dtau_p - dtau_m)/(2*dT) - dtau2)) < 1e-6_pr*maxval(abs(dtau2)), &
      "nrtl_hv_tdep: d2tau/dT2 vs finite differences")

   call nrtl_hv_tdep_linear(T, 0*g, g/R, tau_l, dtau_l, dtau2_l)
   call assert(maxval(abs(tau - tau_l)) < 1e-13_pr*maxval(abs(tau)), "nrtl_hv_tdep: same tau as the linear form A=0, B=g/R")
   call assert(maxval(abs(dtau - dtau_l)) < 1e-13_pr*maxval(abs(dtau)), "nrtl_hv_tdep: same dtau/dT as the linear form")
   call assert(maxval(abs(dtau2 - dtau2_l)) < 1e-13_pr*maxval(abs(dtau2)), "nrtl_hv_tdep: same d2tau/dT2 as the linear form")
end program test_ge_excess_cp
