module yaeos__equilibria_binaries
   !! Module with routines particular to binary mixtures.
   use yaeos__constants, only: pr
   use yaeos__models, only: ArModel

   private
   public :: find_llcl

contains

   subroutine find_llcl(model, z0, zi, P, a, V, T)
      !! # `find_llcl`
      !! Find an initial guess for the critical L-L line of a binary mixture.
      !!
      !! # Description
      !!
      !!
      !! # Examples
      !!
      !!
      !! # References
      !! [1] M. Cismondi, M.L. Michelsen, Global phase equilibrium
      !! calculations:
      !! Critical lines, critical end points and liquid–liquid–vapour
      !! equilibrium in binary mixtures, The Journal of Supercritical Fluids 39
      !!  (2007) 287–295. https://doi.org/10.1016/j.supflu.2006.03.011.
      implicit none
      class(ArModel), intent(in) :: model
      real(kind=pr), intent(in) :: P !! Pressure [bar]
      real(kind=pr), intent(in) :: z0(2) !! Mole fractions of original fluid
      real(kind=pr), intent(in) :: zi(2) !! Mole fractions of new fluid
      real(kind=pr), intent(out) :: a !! Mole fraction of new fluid
      real(kind=pr), intent(out) :: V !! Volume [L/mol]
      real(kind=pr), intent(in out) :: T !! Temperature [K]

      real(kind=pr) :: z(2), lnphi(2), dlnphidn(2, 2)

      integer :: i
      integer :: tries
      real(kind=pr) :: as(50)
      real(kind=pr) :: lambdas(50)

      do tries = 1, 2
         do i = 1, 50
            a = real(i, pr) / 51.
            z = a * zi + (1 - a) * z0
            call model%lnphi_pt(&
                  z, P=P, T=T, V=V, lnPhi=lnPhi, dlnphidn=dlnphidn, root_type=&
                  "liquid")
            lambdas(i) = 1 - dlnphidn(1, 2)

            as(i) = a
         end do

         i = minloc(lambdas, dim=1)
         a = as(i)
         z = a * zi + (1 - a) * z0

         if (lambdas(i) > 0) then
            do while (lambdas(i) > 0)
               T = T - 1
               call model%lnphi_pt(z, P=P, T=T, V=V, lnPhi=lnPhi, dlnphidn=&
                     dlnphidn, root_type="liquid")
               lambdas(i) = 1 - dlnphidn(1, 2)
            end do
         else
            do while (lambdas(i) > 0)
               T = T + 1
               call model%lnphi_pt(z, P=P, T=T, V=V, lnPhi=lnPhi, dlnphidn=&
                     dlnphidn, root_type="liquid")
               lambdas(i) = 1 - dlnphidn(1, 2)
            end do
         end if
      end do
   end subroutine find_llcl
end module yaeos__equilibria_binaries
