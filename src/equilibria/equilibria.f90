module yaeos__equilibria
   !! # Phase Equilibrium Calculations Module
   !!
   !! This module provides comprehensive functionality for calculating
   !! phase equilibria using thermodynamic models. It includes:
   !!
   !! ## Main Features
   !!
   !! - **Flash calculations**: PT, TV, PV flash for vapor-liquid equilibria
   !! - **Saturation points**: Bubble/dew point calculations for mixtures
   !! - **Phase envelopes**: PT, PX, TX envelope tracing
   !! - **Critical points**: Critical point detection and critical lines
   !! - **Stability analysis**: Thermodynamic stability testing
   !! - **Multiphase flash**: Handling systems with more than two phases
   !!
   !! ## Usage Examples
   !!
   !! ### Basic PT Flash
   !! ```fortran
   !! use yaeos, only: flash, EquilibriumState
   !! 
   !! type(EquilibriumState) :: result
   !! call flash(model, n, T=350.0_pr, P=15.0_pr, equilibrium=result)
   !!
   !! if (result%phases == 2) then
   !!     print *, "Vapor fraction:", result%beta(2)
   !! end if
   !! ```
   !!
   !! ### Saturation Pressure
   !! ```fortran
   !! use yaeos, only: saturation_pressure
   !!
   !! real(pr) :: Psat
   !! call saturation_pressure(model, T=350.0_pr, Psat)
   !! ```
   !!
   !! ## Available Procedures
   !!
   !! - [[flash]]: Main flash calculation routine
   !! - [[saturation_pressure]]: Bubble/dew pressure calculations  
   !! - [[saturation_temperature]]: Bubble/dew temperature calculations
   !! - [[pt_envelope_2ph]]: PT phase envelope tracing
   !! - [[px_envelope_2ph]]: PX phase envelope tracing
   !! - [[critical_point]]: Critical point calculation
   !! - [[min_tpd]]: Stability analysis using TPD method
   !!
   !! ## Derived Types
   !!
   !! - [[EquilibriumState]]: Contains results of equilibrium calculations
   !! - [[MPEquilibriumState]]: Results for multiphase equilibria
   !! - [[PTEnvelope2Ph]]: PT envelope data structure
   !! - [[CriticalLine]]: Critical line data structure
   !!
   !! ## References
   !!
   !! 1. Michelsen, M. L., & Mollerup, J. M. (2007). Thermodynamic models:
   !!    Fundamentals & computational aspects. Tie-Line Publications.
   !! 2. Firoozabadi, A. (1999). Thermodynamics of hydrocarbon reservoirs.
   !!    McGraw-Hill Professional.

   ! Stability analysis
   use yaeos__equilibria_stability, only: tm, min_tpd

   ! Equilibrium State definitions
   use yaeos__equilibria_equilibrium_state, only: &
      EquilibriumState, MPEquilibriumState

   ! Phase split calculations
   use yaeos__equilibria_flash, only: flash

   use yaeos__equilibria_multiphase_flash, only: &
      solve_mp_flash_point, pt_mp_flash, MPEquilibriumState

   ! Saturation points
   use yaeos__equilibria_saturation_points, only:&
      saturation_pressure, saturation_temperature

   ! Critical points
   use yaeos__equilibria_critical, only: &
      critical_line, CriticalLine, critical_point, spec_CP

   use yaeos__equilibria_binaries, only: &
      find_llcl

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
   use yaeos__equilibria_boundaries_generalized_isopleths, only: create_generalized_isoz_line, GeneralizedIsoZLine
end module yaeos__equilibria
