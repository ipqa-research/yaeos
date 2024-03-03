module yaeos_phase_equilibria_stability
   use yaeos_constants, only: pr
   use yaeos_models_ar, only: ArModel
   use yaeos_thermoprops, only: fugacity_vt
   implicit none

contains

   subroutine tangent_plane_distance(z, lnphi_z, w, lnphi_w, tpd)
      real(pr), intent(in) :: z(:) !! Test phase composition
      real(pr), intent(in) :: lnphi_z(:) !! Test phase composition
      real(pr), intent(in) :: w(:) !! Trial phase composition
      real(pr), intent(in) :: lnphi_w(:) !! Trial phase composition

      real(pr), intent(out) :: tpd !! Tangent plane distance

      tpd = sum(w * (log(w) + lnphi_w - log(z) - lnphi_z))
   end subroutine
end module
