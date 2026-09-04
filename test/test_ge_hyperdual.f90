module nrtl_hd
   !! Example implementation of the NRTL model using automatic differentiation.
   !!
   !! The NRTL model is represented by the equations:
   !!
   !! \[tau = A + B/T\]
   !! \[G = exp(-C * tau)\]
   !! \[Ge = R * T * sum(n * sum(n * tau * G) / sum(n * G))\]
   use yaeos
   use yaeos__adiff_hyperdual_ge_api, only: GeModelAdiff
   use hyperdual_mod
   implicit none

   type, extends(GeModelAdiff) :: NRTLHD
      !! NRTL model with automatic differentiation.
      real(pr), allocatable :: A(:,:), B(:,:), C(:,:)
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
      G = exp(-self%C * tau)

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

   real(pr) :: n(3), T

   real(pr) :: Geog, Geimp
   real(pr) :: GeogT, GeimpT
   real(pr) :: GeogT2, GeimpT2
   real(pr) :: Geogn(3), Geimpn(3)
   real(pr) :: GeognT(3), GeimpNt(3)
   real(pr) :: Geogn2(3,3), Geimpn2(3,3)

   real(pr) :: A(3,3), B(3,3), C(3,3)

   A = 0; B=0; C=0   

   A(1, 2) = 0.1; A(1, 3) = 0.2
   A(2, 1) = 0.3; A(2, 3) = 0.4
   A(3, 1) = 0.5; A(3, 2) = 0.6

   B(1, 2) = 300; B(1, 3) = 500
   B(2, 1) = 700; B(2, 3) = 900
   B(3, 1) = 1100; B(3, 2) = 1300

   C(1, 2) = 0.7; C(1, 3) = 0.8; C(2, 3) = 1.0
   C(2, 1) = 1.1; C(3, 1) = 1.2; C(3, 2) = 1.3

   og%a = A
   og%b = B
   og%C = C

   imp%A = og%a
   imp%B = og%b
   imp%C = og%C

   n = [0.2, 0.8, 0.1]
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