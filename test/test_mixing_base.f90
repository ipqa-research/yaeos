program main
   use yaeos__models_ar_cubic_mixing_base
   implicit none

   integer, parameter :: nc=3
   real(pr) :: d1i(nc), d1, dd1i(nc), dd1ij(nc, nc)
   real(pr) :: L, dLi(nc), dLij(nc, nc)
   real(pr) :: dLi_num(nc), dlij_num(nc, nc)
   real(pr) :: n(nc), dn, f1, f2, f3, f4

   integer :: i, j

   n = [1, 5, 2]
   d1i = [0.1_pr, 0.2_pr, 0.3_pr]
   call d1mix_rkpr(n, d1i, d1, dd1i, dd1ij)
   call lamdba_hv(d1, dd1i, dd1ij, L, dLi, dLij)

   print *, "L  = ", L

   print *, 'd1i  = ', d1i
   print *, 'dd1i = ', dd1i

   dn = 0.0001_pr
   do i=1,nc
      n(i) = n(i) + dn
      dLi_num(i) = (Lval(n) - L)/dn
      n(i) = n(i) - dn
      do j=1, nc
         n(i) = n(i) + dn
         n(j) = n(j) + dn 
         f1 = Lval(n)
         n(j) = n(j) - 2*dn
         f2 = Lval(n)
         n(i) = n(i) - 2*dn
         f3 = Lval(n)
         n(j) = n(j) + 2*dn
         f4 = Lval(n)
         dlij_num(i, j) = (f1 - f2 + f3 - f4)/(4*dn**2)
      end do
   end do

   call d1mix_rkpr(n, d1i, d1, dd1i, dd1ij)
   call lamdba_hv(d1, dd1i, dd1ij, L, dLi, dLij)

   print *, 'dLi_num = ', dLi_num
   print *, 'dLi     = ', dLi

   print *, 'dlij_num = ', dlij_num
   print *, 'dlij     = ', dLij

contains

   real(pr) function Lval(n)
      real(pr), intent(in) :: n(:)

      call d1mix_rkpr(n, d1i, d1, dd1i, dd1ij)
      call lamdba_hv(d1, dd1i, dd1ij, Lval, dLi, dLij)
   end function

end program