module yaeos__equilibria
   use yaeos__equilibria_equilibria_state, only: EquilibriumState
   use yaeos__equilibria_flash, only: flash
   use yaeos__equilibria_saturation_points, only:&
       saturation_pressure, saturation_temperature
   use yaeos__phase_equilibria_boundaries_phase_envelopes_pt, only:&
       PTEnvel2, pt_envelope_2ph
   use yaeos__phase_equilibria_boundaries_phase_envelopes_px, only:&
       PXEnvel2, px_envelope_2ph
   use yaeos__phase_equilibria_auxiliar, only: k_wilson
   implicit none
end module