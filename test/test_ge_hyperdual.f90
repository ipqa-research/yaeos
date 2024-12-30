module nrtl_hd
   !! Example implementation of the NRTL model using automatic differentiation.
   !!
   !! The NRTL model is represented by the equations:
   !!
   !! \[tau = A + B/T\]
   !! \[G = exp(-alpha * tau)\]
   !! \[Ge = R * T * sum(n * sum(n * tau * G) / sum(n * G))\]
   use yaeos
   use yaeos__adiff_hyperdual_ge_api, only: GeModelAdiff
   use hyperdual_mod
   implicit none

   type, extends(GeModelAdiff) :: NRTLHD
      !! NRTL model with automatic differentiation.
      real(pr), allocatable :: A(:,:), B(:,:)
      real(pr) :: alpha
   contains
      procedure :: Ge => Ge
   end type NRTLHD

contains
   function ge(self, n, t)
      class(NRTLHD) :: self
      type(hyperdual), intent(in) :: n(:)
      type(hyperdual), intent(in) :: t
      type(hyperdual) :: ge

      type(hyperdual) ::tau(size(n), size(n)), G(size(n), size(n))
      

      real(pr) :: down
      integer :: i, j

      tau = self%a(:, :) + self%b(:, :)/T
      G = exp(-self%alpha * tau)

      ge = 0._pr
      do i=1,size(n)
         ge = ge + n(i) * sum(n(:) * tau(:, i) * G(:, i)) / sum(n(:) *  g(:, i))
      end do
      ge = R * T * ge
   end function ge
end module nrtl_hd

program test
   use yaeos
   use nrtl_hd
   use fixtures_models, only: binary_NRTL_tape
   use testing_aux, only: test_title, assert

   type(NRTL) :: og
   type(NRTLHD) :: imp

   real(pr) :: n(2), T

   real(pr) :: Geog, Geimp
   real(pr) :: GeogT, GeimpT
   real(pr) :: GeogT2, GeimpT2
   real(pr) :: Geogn(2), Geimpn(2)
   real(pr) :: GeognT(2), GeimpNt(2)
   real(pr) :: Geogn2(2,2), Geimpn2(2,2)

   og = binary_NRTL_tape()

   imp%A = og%a
   imp%B = og%b
   imp%alpha = og%C(1,2)

   n = [0.2, 0.8]
   T = 273

   write(*, *) test_title("NRTL model with automatic differentiation")

   call og%excess_gibbs(n, T, Geog)
   call imp%excess_gibbs(n, T, Geimp)

   call assert(abs(Geog - Geimp) < 1e-5, "Ge")

   call og%excess_gibbs(n, T, Geog, GeT=GeogT)
   call imp%excess_gibbs(n, T, Geimp, GeT=GeimpT)

   call assert(abs(GeogT - GeimpT) < 1e-5, "GeT")

   call og%excess_gibbs(n, T, Ge=Geog, GeT2=GeogT2)
   call imp%excess_gibbs(n, T, Ge=Geimp, GeT2=GeimpT2)

   call assert(abs(GeogT2 - GeimpT2) < 1e-5, "GeT2")

   call og%excess_gibbs(n, T, Ge=Geog, GeTn=GeognT)
   call imp%excess_gibbs(n, T, Ge=Geimp, GeTn=GeimpnT)

   call assert(all(abs(GeognT - GeimpnT) < 1e-5), "GeTn")

   call og%excess_gibbs(n, T, Geog, Gen=Geogn)
   call imp%excess_gibbs(n, T, Geimp, Gen=Geimpn)

   call assert(all(abs(Geogn - Geimpn) < 1e-5), "Gen")

   call og%excess_gibbs(n, T, Geog, Gen2=Geogn2)
   call imp%excess_gibbs(n, T, Geog, Gen2=Geimpn2)

   call assert(all(abs(Geogn2 - Geimpn2) < 1e-5), "Gen2")

end program