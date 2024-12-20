module yaeos__equilibria
   !! Module to handle phase equilibria calculations.

   ! Stability analysis
   use yaeos__equilibria_stability, only: tm, min_tpd

   ! Equilibrium State definitions
   use yaeos__equilibria_equilibrium_state, only: EquilibriumState

   ! Pure component saturation pressure
   use yaeos__equilibria_pure_psat, only: Psat

   ! Phase split calculations
   use yaeos__equilibria_flash, only: flash

   ! Saturation points
   use yaeos__equilibria_saturation_points, only:&
      saturation_pressure, saturation_temperature

   ! Critical points
   use yaeos__equilibria_critical, only: &
      critical_line, CriticalLine, critical_point, spec_CP

   ! Phase equilibria boundaries
   use yaeos__equilibria_boundaries_phase_envelopes_pt, only:&
      PTEnvel2, pt_envelope_2ph, find_hpl
   use yaeos__equilibria_boundaries_phase_envelopes_px, only:&
      PXEnvel2, px_envelope_2ph

   ! Extra
   use yaeos__equilibria_auxiliar, only: k_wilson, p_wilson
   implicit none
end module yaeos__equilibria
