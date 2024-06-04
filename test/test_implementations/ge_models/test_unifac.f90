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
   real(pr) :: theta(ng), dthetadx(ng, nc)
   real(pr) :: lngamma(nc), dlngamma_dn(nc, nc)
   real(pr) :: ln_Gamma(ng), dln_Gammadt(ng)=0, dln_Gammadt_num(ng)=0, dln_Gammadn(ng, nc)


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

   call group_big_gamma(model, x, T, ln_Gamma=ln_Gamma)
   print *, "ln_Gamma: ", maxval(abs(ln_Gamma - ln_Gamma_val)) < 1e-7
   
   call model%ln_activity_coefficient(x, T, lngamma)
   print *, "lngamma: ", maxval(abs(lngamma - lngamma_val)) < 1e-7

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
