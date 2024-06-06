program main
   use iso_fortran_env, only: int64
   use yaeos, only: pr
   use yaeos__models_ge_group_contribution_unifac
   use stdlib_io_npy, only: load_npy
   implicit none

   integer :: i

   type(UNIFAC) :: model
   integer, parameter :: nc = 3, ng = 4
   real(pr) :: x(nc) = [0.2, 0.7, 0.1], T=150

   ! integer, parameter :: nc = 2, ng = 3
   ! real(pr) :: x(nc) = [0.3, 0.7], T=150

   real(pr), allocatable :: Aij(:, :)
   real(pr), allocatable :: Qk(:), Rk(:)
   real(pr) :: dx=1e-5, dpsidt_num(ng, ng), dpsidt(ng, ng)

   type(Groups) :: molecules(nc)
   real(pr) :: psi(ng, ng)
   real(pr) :: theta(ng), theta_ji(ng, nc), dthetadx(ng, nc)
   real(pr) :: lngamma(nc), dlngamma_dn(nc, nc)
   real(pr) :: ln_Gamma(ng), dln_Gammadt(ng)=0, dln_Gammadt_num(ng)=0, dln_Gammadn(ng, nc)
   real(pr) :: ln_gamma_c(nc), dln_gamma_c_dn(nc, nc)

   real(pr) :: Ge, Gen(nc)


   real(pr) :: lngamma_val(nc) = [0.84433780935070013, -0.19063836609197171, -2.9392550019369406]
   real(pr) :: ln_Gamma_val(ng) = [0.43090864639734738, 0.27439937388510327,  0.52442445057961795, -2.8793657040300329]

   integer(int64) :: rate, st, et
   call system_clock(count_rate=rate)

   call load_npy("data/unifac_aij.npy", Aij)
   call load_npy("data/unifac_Qk.npy", Qk)
   call load_npy("data/unifac_Rk.npy", Rk)

   ! ! Ethane [CH3]
   molecules(1)%groups_ids = [1]
   molecules(1)%number_of_groups = [2]

   ! ! Ethanol [CH3, CH2, OH]
   molecules(2)%groups_ids = [1, 2, 14]
   molecules(2)%number_of_groups = [1, 1, 1]

   ! ! Methylamine [H3C-NH2]
   molecules(3)%groups_ids = [28]
   molecules(3)%number_of_groups = [1]


   model = setup_unifac(&
      molecules, &
      Eij=Aij, &
      Qk=Qk, &
      Rk=Rk &
      )

   print *, "=================================================================="
   print *,  "Cosas que dan"
   print *, "=================================================================="
   ! Coeficientes de actividad
   call model%ln_activity_coefficient(x, T, lngamma)
   print *, "lngamma: ", maxval(abs(lngamma - lngamma_val)) < 1e-7
   print *, " "

   ! Valor de los big gamma del residual
   call group_big_gamma(model, x, T, ln_Gamma=ln_Gamma)
   print *, "ln_Gamma: ", maxval(abs(ln_Gamma - ln_Gamma_val)) < 1e-7
   print *, " "

   ! Ge
   call model%excess_gibbs(x, T, Ge=Ge, Gen=Gen)
   print *, "Ge: ", Ge
   print *, " "

   ! Gen
   print *, "Gen: ", Gen
   print *, " "

   ! Derivada composicional del gamma composicional
   call combinatorial_activity(model, x, ln_gamma_c=ln_gamma_c, dln_gamma_c_dn=dln_gamma_c_dn)
   x = [20, 70, 10]
   print *, x
   print *, "ln_gamma_c: ", ln_gamma_c, "|| expected:", [-0.017479483844929714, -0.004513178803947096, -0.059889298605798724]
   print *, " "
   print *, "dln_gamma_c_dn: "
   print *,  dln_gamma_c_dn
   print *, "expected:"
   print * , &
      [0.031692694440355496, -0.01793757017123132, 0.06217760231790814, &
      -0.017937570171231723, 0.009591847252321026, -0.03126779042378376, &
      0.0621776023179087, -0.03126779042378558, 0.09451932833068105]
   print *, " "

   ! Thetas totales
   call thetas(model, x, theta_j=theta)
   print *, "Total Thetas: ", theta
   print *, "Expected: ", [0.40465035571750835, 0.16397709526288395, 0.3643935450286309, 0.06697900399097695]
   print *, " "

   ! Thetas_i
   ! call thetas_i(model, theta_ji)
   ! print *, "Thetas_i: "
   ! print *, theta_ji(1,1), theta_ji(1,2), theta_ji(1,3)
   ! print *, theta_ji(2,1), theta_ji(2,2), theta_ji(2,3)
   ! print *, theta_ji(3,1), theta_ji(3,2), theta_ji(3,3)
   ! print *, theta_ji(4,1), theta_ji(4,2), theta_ji(4,3)
   ! print *, "Expected: "
   ! print *, [1., 0.32766615, 0.]
   ! print *, [0., 0.20865533, 0.]
   ! print *, [0., 0.46367852, 0.]
   ! print *, [0., 0.        , 1.]
   ! print *, " "

   ! ==========================================================================
   ! Cosas que no dan
   ! --------------------------------------------------------------------------
   call group_big_gamma(model, x, T, ln_Gamma, dln_gammadt=dln_Gammadt)
   call group_big_gamma(model, x, T+dx, dln_Gammadt_num)

   print *, "dlnGamma_dT:"
   print *, "numm: ", (dln_Gammadt_num - ln_gamma)/dx
   print *, "anal: ", dln_gammadt

   print *, "dlnGamma_dx"
   call group_big_gamma(model, x, T, ln_gamma, dln_Gammadx=dln_Gammadn)
   print *, dln_Gammadn

   call residual_activity(model, x, T, lngamma, dln_gamma_dn=dlngamma_dn)

end program main
