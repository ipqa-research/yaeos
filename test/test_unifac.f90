program main
   use yaeos, only: pr
   use yaeos_models_ge_group_contribution_unifac
   use stdlib_io_npy, only: load_npy
   implicit none

   integer, parameter :: nc = 3, ng = 4
   integer :: i

   type(Groups) :: molecules(nc)
   type(UNIFAC) :: model
   real(pr) :: x(nc) = [0.2, 0.7, 0.1]

   real(pr), allocatable :: Aij(:, :)
   real(pr), allocatable :: Qk(:), Rk(:)

   real(pr) :: psi(ng, ng)
   real(pr) :: theta(ng), dthetadx(ng, nc)
   real(pr) :: lngamma(nc)

   call load_npy("data/unifac_aij.npy", Aij)
   call load_npy("data/unifac_Qk.npy", Qk)
   call load_npy("data/unifac_Rk.npy", Rk)

   ! Ethane [CH3]
   molecules(1)%groups_ids = [1]
   molecules(1)%number_of_groups = [2]
   
   ! Ethanol [CH3, CH2, OH]
   molecules(2)%groups_ids = [1, 2, 14]
   molecules(2)%number_of_groups = [1, 1, 1]
   
   ! Methylamine [H3C-NH2]
   molecules(3)%groups_ids = [28]
   molecules(3)%number_of_groups = [1]

   model = setup_unifac(&
      molecules, &
      Eij=Aij, &
      Qk=Qk, &
      Rk=Rk &
   )

   psi = 0
   call model%psi_function%psi(model%groups_stew, 200._pr, psi)
   call excess_gibbs(model, x, 200._pr)

   call group_area_fraction(model, x, theta, dthetadx)
   print *, theta

   print *, "Derivs:"
   do i=1,ng
      print *, dthetadx(i, :)
   end do

   call ln_activity_coefficient(model, x, 200.0_pr, lngamma)
   print *, "lngamma: "
   print *, lngamma
   !  0.80338267603153490       -9.7236411677591672E-002 -0.30836537797129104     

   call group_big_gamma(model, X, 200._pr, theta, dthetadx)
   print *, ""
   print *, theta
   print *, ""
   do i=1,ng
      print *, dthetadx(i, :)
   end do
end program main
