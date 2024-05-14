program main
   use iso_fortran_env, only: int64
   use yaeos, only: pr
   use yaeos_models_ge_group_contribution_unifac
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
   real(pr) :: lngamma(nc)
   real(pr) :: ln_Gamma(ng), dln_Gammadt(ng)=0, dln_Gammadt_num(ng)=0

   integer(int64) :: rate, st, et
   call system_clock(count_rate=rate)

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

   ! call model%psi_function%psi(model%groups_stew, T, psi, dpsidt=dpsidt)
   ! print *, "psi"
   ! do i=1,ng
   !    print *, dpsidt(:, i)
   ! end do

   block
      integer :: i, j
      associate (ids => model%groups_stew%groups_ids)
      do i=1,size(model%groups_stew%groups_ids)
         exit
         do j=1,size(model%groups_stew%groups_ids)
            print *, ids(i),ids(j), Aij(ids(i), ids(j))
         end do
      end do
      end associate
   end block


   call model%ln_activity_coefficient(x, T, lngamma)
   print *, "lngamma: ", lngamma
   call ln_activity_coefficient(model, x, T, lngamma)
   print *, "lngamma: ", lngamma

   call group_big_gamma(model, x, T, ln_Gamma, dln_gammadt=dln_Gammadt)
   call group_big_gamma(model, x, T+dx, dln_Gammadt_num)

   print *, "numm: ", (dln_Gammadt_num - ln_gamma)/dx
   print *, "anal: ", dln_gammadt

   call group_big_gamma(model, x, T, ln_Gamma=ln_Gamma)
   print *, "ln_Gamma: ", ln_Gamma
   do i=1, nc
      x = 0
      x(i) = 1
      call group_big_gamma(model, x, T, ln_Gamma)
      print *, "ln_gamma_pure:", ln_gamma
   end do
end program main
