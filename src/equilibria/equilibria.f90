module yaeos__equilibria
   !! Module to handle phase equilibria calculations.

   ! Stability analysis
   use yaeos__equilibria_stability, only: tm, min_tpd

   ! Equilibrium State definitions
   use yaeos__equilibria_equilibrium_state, only: EquilibriumState

   ! Phase split calculations
   use yaeos__equilibria_flash, only: flash

   ! Saturation points
   use yaeos__equilibria_saturation_points, only:&
      saturation_pressure, saturation_temperature

   ! Critical points
   use yaeos__equilibria_critical, only: &
      critical_line, CriticalLine, critical_point, spec_CP

   ! Extra
   use yaeos__equilibria_auxiliar, only: k_wilson, p_wilson

   ! Phase equilibria boundaries
   use yaeos__equilibria_boundaries_pure_saturation, only: &
      PurePsat, pure_saturation_line
   use yaeos__equilibria_boundaries_phase_envelopes_px, only: &
      PXEnvel2, px_envelope_2ph
   use yaeos__equilibria_boundaries_phase_envelopes_tx, only: &
      TXEnvel2, tx_envelope_2ph
   use yaeos__equilibria_boundaries_phase_envelopes_pt, only: &
      PTEnvel2, pt_envelope_2ph, find_hpl
   use yaeos__equilibria_boundaries_phase_envelopes_pt3, only: &
      PTEnvel3, pt_envelope_3ph
   use yaeos__equilibria_boundaries_phase_envelopes_px3, only: &
      PXEnvel3, PX_envelope_3ph
   use yaeos__equilibria_boundaries_phase_envelopes_mp, only: &
      PTEnvelMP, pt_envelope
   use yaeos__equilibria_boundaries_phase_envelopes_mp_px, only: &
      PXEnvelMP, px_envelope
end module yaeos__equilibria
