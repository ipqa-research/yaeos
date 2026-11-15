program test_flash_ge
   use yaeos
   use testing_aux, only: assert, test_title
   implicit none

   integer, parameter :: nc = 3
   type(UNIFAC) :: model
   type(Groups) :: molecules(nc)
   real(pr) :: n(nc), w(nc), T
   type(EquilibriumState) :: fr
   real(pr) :: mintpd
   integer :: i


   write(*, *) test_title("Testing FLASH GE")
   
   T = 250
   n = [0.2, 0.7, 0.1]

   molecules(1)%groups_ids = [1, 42]
   molecules(1)%number_of_groups = [1, 1]

   molecules(2)%groups_ids = [1, 2]
   molecules(2)%number_of_groups = [2, 4]

   molecules(3)%groups_ids = [16]
   molecules(3)%number_of_groups = [1]

   ! setup UNIFAC model
   model = setup_unifac(molecules)

   call min_tpd(model, n, P=1._pr, T=T, mintpd=mintpd, w=w)
   
   call assert(mintpd < 0.0_pr, "Reach a negative tm")

   fr = flash(model, n, T, k0=w/n, iters=i)

   call assert(fr%kind /= "failed", "Should converge")
   call assert(abs(fr%beta - 0.2636) < 1e-4_pr, "beta value")
   call assert(maxval(abs(fr%x - [0.0568, 0.9413, 0.0018])) < 1e-4_pr, "x phase")
   call assert(maxval(abs(fr%y - [0.5999, 0.02569,  0.3743])) < 1e-4_pr, "y phase")
end program test_flash_ge
