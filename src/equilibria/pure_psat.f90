module yaeos__equilibria_pure_psat
   !! Module used to calculate the saturation pressure of pure components at
   !! a given temperature.
   use yaeos__constants, only: pr
   use yaeos__models, only: ArModel, size
contains
   real(pr) function Psat(eos, ncomp, T)
      !! Calculation of saturation pressure of a pure component using the
      !! secant method.
      class(ArModel), intent(in) :: eos !! Model that will be used
      integer, intent(in) :: ncomp 
         !! Number of component in the mixture from which the saturation pressure
         !! will be calculated
      real(pr), intent(in) :: T !! Temperature [K]

      real(pr) :: P1, P2
      real(pr) :: f1, f2

      real(pr) :: n(size(eos))

      n = 0
      n(ncomp) = 1

      P1 = 0.5
      P2 = 1

      do while(abs(diff(P2)) > 1e-5)
         f1 = diff(P1)
         f2 = diff(P2)
         Psat = (P1 * f2 - P2 * f1)/(f2 - f1)
         P1 = P2
         P2 = Psat
      end do
   contains
      real(pr) function diff(P)
         real(pr), intent(in) :: P
         real(pr) :: V_l, V_v
         real(pr) :: phi_v(size(eos)), phi_l(size(eos))
         call eos%lnphi_pt(n, P, T, V=V_v, lnPhi=phi_v, root_type="vapor")
         call eos%lnphi_pt(n, P, T, V=V_l, lnPhi=phi_l, root_type="liquid")
         diff = phi_v(ncomp) - phi_l(ncomp)
      end function
   end function Psat
end module