module yaeos__equilibria
   !! # Phase-equilibrium API
   !! Public entry points for stability, saturation, flash, critical-point, and
   !! binary three-phase calculations.
   !!
   !! # Description
   !! This façade re-exports the principal equilibrium result types and
   !! calculation routines from their implementation modules. Importing
   !! `yaeos` is usually more convenient; this module is useful when a client
   !! wants only the equilibrium API.
   !!
   !! # Examples
   !! ```fortran
   !! use yaeos__equilibria, only: EquilibriumState, flash
   !! result = flash(model, z, T, p_spec=P)
   !! ```

   ! Stability analysis
   use yaeos__equilibria_stability, only: tm, min_tpd

   ! Equilibrium State definitions
   use yaeos__equilibria_equilibrium_state, only: &
      EquilibriumState, MPEquilibriumState

   ! Phase split calculations
   use yaeos__equilibria_flash, only: flash, flash_no_beta_limits
   use yaeos__equilibria_flash_second_order, only: flash_second_order

   use yaeos__equilibria_multiphase_flash, only: &
      solve_mp_flash_point, pt_mp_flash, MPEquilibriumState

   ! Saturation points
   use yaeos__equilibria_saturation_points, only:&
      saturation_pressure, saturation_temperature

   ! Critical points
   use yaeos__equilibria_critical, only: &
      critical_line, CriticalLine, critical_point, spec_CP

   use yaeos__equilibria_binaries, only: &
      find_llcl, BinaryThreePhase, binary_llv_from_cep

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
   use yaeos__equilibria_boundaries_phase_envelopes_mp_tx, only: &
      TXEnvelMP, tx_envelope
   use yaeos__equilibria_boundaries_generalized_isopleths, only: &
      create_generalized_isoz_line, GeneralizedIsoZLine
end module yaeos__equilibria
