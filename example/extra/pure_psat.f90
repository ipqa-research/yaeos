!! Program to calculate the vapor pressure of pure components
!!
program main
   use yaeos
   use forsus, only: Substance, forsus_default_dir, forsus_dir
   implicit none
   integer, parameter :: nc=2
   type(CubicEoS) :: eos
   type(Substance) :: sus(nc)
   real(pr) :: n(nc), T, x, f, df
   integer :: i

   forsus_dir = "build/dependencies/forsus/" // forsus_default_dir
   sus(1) = Substance("water")
   sus(2) = Substance("ethanol")

   eos = SoaveRedlichKwong(&
      sus%critical%critical_temperature%value, &
      sus%critical%critical_pressure%value/1e5,&
      sus%critical%acentric_factor%value &
      )

   n = [1, 0]

   T = 273.15_pr + 50

   do i=273+90, nint(sus(1)%critical%critical_temperature%value)
      T = real(i,pr)
      call foo(T, f)
      print *, T, F
   end do

contains

   subroutine foo(x, f)
      real(pr), intent(in) :: x
      real(pr), intent(out) :: f
      real(pr) :: phi_v(nc), phi_l(nc), dphidV_v, dphidV_l, Ar, P
      real(pr) :: dlnphi_dt_v(nc), dlnphi_dt_l(nc)
      real(pr) :: V_l, V_v
      real(pr) :: diff

      P = 1
      diff = 1
      do while(abs(diff) > 1e-8_pr .and. .not. isnan(diff))
         call eos%lnphi_pt(n, P, x, V=V_v, lnPhi=phi_v, dlnPhidt=dlnphi_dt_v, root_type="vapor")
         call eos%lnphi_pt(n, P, x, V=V_l, lnPhi=phi_l, dlnPhidt=dlnphi_dt_l, root_type="liquid")
         diff = exp(phi_v(1)) - exp(phi_l(1))
         P = P + sign(diff, 0.1_pr)
      end do

      f = P
   end subroutine
end program main
