program main
   !! Test the calculation of fugacity. Compared to numerical derivatives of 
   !! manual calculations using the fugacity coefficient routine.
   use yaeos
   use fixtures_models, only: binary_PR76
   implicit none
   class(ArModel), allocatable :: model

   integer, parameter :: nc = 2

   real(pr) :: n(nc), V, T, P, lnf(nc), dlnfdt(nc), dlnfdv(nc), dlnfdn(nc,nc)
   real(pr) :: numdiff(nc), numdiff_dn(nc, nc), dx, dn(nc), lnfdn(nc)

   integer :: i
   
   model = binary_PR76()
   n = [8, 5]
   T = 200._pr
   V = model%get_v0(n, P, T)*5

   call model%lnfug_vt(&
      n=n, V=V, T=T, lnf=lnf, &
      P=P, dlnfdt=dlnfdt, &
      dlnfdv=dlnfdv, dlnfdn=dlnfdn)

   print *, "lnf"
   print *, lnf
   print *, lnfug(n, V, T)

   dx = 0.0001
   print *, "dt"
   numdiff = (lnfug(n, V, T+dx) - lnfug(n, V, T))/dx
   print *, dlnfdt
   print *, numdiff
   
   print *, "dv"
   dx = 0.0001_pr
   numdiff = (lnfug(n, V+dx, T ) - lnfug(n, V, T))/dx
   print *, dlnfdv
   print *, numdiff

   print *, "dn"
   do i=1,nc
      dn = 0.0_pr
      dn(i) = 0.000001_pr
      numdiff_dn(i, :) = (lnfug(n+dn, V, T) - lnfug(n, V, T))/dn(i)
   end do

   print *, dlnfdn
   print *, numdiff_dn

contains

   function lnfug(n, V, T)
      real(pr), intent(in) :: n(:), V, T
      real(pr) :: lnfug(nc)

      real(pr) :: lnphi(nc), P

      call model%lnphi_vt(n, V, T, lnphi=lnphi, P=P)
      lnfug = log(n/sum(n)) + lnphi + log(P)
   end function
end program main