!! Program to calculate the vapor pressure of pure components
module pure_psat
   !! Module used to calculate the saturation pressure of pure components at
   !! a given temperature.
   use yaeos, only: ArModel, pr, R, size
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

      do while(abs(diff(P2)) > 1e-10)
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
         diff = phi_v(1) - phi_l(1)
      end function
   end function Psat
end module

program main
   use yaeos, only: CubicEoS, SoaveRedlichKwong, pr
   use forsus, only: Substance, forsus_default_dir, forsus_dir
   use pure_psat, only: Psat
   
   implicit none

   integer, parameter :: nc=2
   type(CubicEoS) :: eos
   type(Substance) :: sus(nc)
   real(pr) :: n(nc), T
   integer :: i, j

   forsus_dir = "build/dependencies/forsus/" // forsus_default_dir
   sus(1) = Substance("water")
   sus(2) = Substance("ethanol")
   eos = SoaveRedlichKwong(&
      sus%critical%critical_temperature%value, &
      sus%critical%critical_pressure%value/1e5,&
      sus%critical%acentric_factor%value &
      )


   T = 273.15_pr + 50
   do i=273+90, nint(maxval(sus%critical%critical_temperature%value))
      T = real(i,pr)
      print *, T, (Psat(eos, j, T), j=1,nc)
   end do
end program main
