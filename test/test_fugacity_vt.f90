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


   if (sum(abs(lnf - lnfug(n, V, T))) > 1e-6) then
      error stop "lnf failed"
   end if

   dx = 0.0001
   numdiff = (lnfug(n, V, T+dx) - lnfug(n, V, T))/dx

   if (sum(abs(dlnfdt - numdiff)) > 1e-6) then
      error stop "dlnfdt failed"
   end if
   
   dx = 0.0001_pr
   numdiff = (lnfug(n, V+dx, T ) - lnfug(n, V, T))/dx
   if (sum(abs(dlnfdv - numdiff)) > 1e-4) then
      error stop "dlnfdv failed"
   end if

   do i=1,nc
      dn = 0.0_pr
      dn(i) = 0.000001_pr
      numdiff_dn(i, :) = (lnfug(n+dn, V, T) - lnfug(n, V, T))/dn(i)
   end do

   if (sum(abs(dlnfdn - numdiff_dn)) > 1e-6) then
      error stop "dlnfdn failed"
   end if

contains

   function lnfug(n, V, T)
      real(pr), intent(in) :: n(:), V, T
      real(pr) :: lnfug(nc)

      real(pr) :: lnphi(nc), P

      call model%lnphi_vt(n, V, T, lnphi=lnphi, P=P)
      lnfug = log(n/sum(n)) + lnphi + log(P)
   end function
end program main